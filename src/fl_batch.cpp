#include "fl_batch.hpp"

#include "fl_field.hpp"
#include "fl_image.hpp"
#include "fl_io.hpp"
#include "fl_project.hpp"
#include "fl_raster.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#include <cstddef>

namespace flatland {

/* ----------------------
   Output Helpers
   ---------------------- */
std::string esc(const std::string& s) {
    std::string out;
    for(char c:s) (c=='\\'||c=='"') ? out+="\\" + std::string(1,c) : out+=c;
    return out;
}

// Emit a JSON-safe number: non-finite values become null so pipelines don't choke.
template <typename T>
std::string jnum(T v) {
    if (!std::isfinite((double)v)) return "null";
    std::ostringstream o; o << std::setprecision(10) << (double)v; return o.str();
}

/* ----------------------
   Application Runner
   ---------------------- */
template <typename T>
void run_app(const std::string& obj, const std::string& default_data, const std::string& out_pre,
             const std::vector<BatchEntry>& batch, double default_res, bool cull, bool json,
             unsigned threads, ValueMode forced_mode) {

    Mesh<T> mesh = load_mesh<T>(obj);

    // Each distinct field file is loaded ONCE into a shared, read-only cache; workers
    // then extract their timestep's column from it concurrently. This keeps the
    // common time-series case (one matrix file shared by every view) from being
    // re-parsed per thread, while distinct per-view files are still each loaded once.
    std::map<std::string, std::shared_ptr<FieldMatrix<T>>> mat_cache;
    std::mutex cache_mtx;
    auto get_matrix = [&](const std::string& path) -> std::shared_ptr<FieldMatrix<T>> {
        std::lock_guard<std::mutex> lk(cache_mtx);
        auto it = mat_cache.find(path);
        if (it == mat_cache.end()) {
            auto m = std::make_shared<FieldMatrix<T>>();
            load_matrix_into(path, mesh, *m, forced_mode);
            it = mat_cache.emplace(path, m).first;
        }
        return it->second;
    };
    // Pre-warm + validate the default field up front so a bad -d fails fast.
    if (!default_data.empty()) get_matrix(parse_field_token(default_data).path);

    std::vector<T> view_res(batch.size());
    for (size_t k=0; k<batch.size(); ++k)
        view_res[k] = (batch[k].resolution > 0) ? (T)batch[k].resolution : (T)default_res;

    std::vector<ViewResult<T>> results(batch.size());
    std::atomic<size_t> next{0};
    std::mutex err_mtx;
    std::string err_msg;

    // Worker: pull timestep indices off a shared counter. Each owns a renderer and a
    // reusable field buffer, so per-timestep field files load on demand -> peak memory
    // is O(threads), not O(timesteps). The shared fixed field is never reloaded.
    auto worker = [&]() {
        Renderer<T> renderer;
        Field<T> local;            // column extracted for the current view
        size_t k;
        while ((k = next.fetch_add(1)) < batch.size()) {
            try {
                const auto& b = batch[k];
                const std::string& token = b.data_file.empty() ? default_data : b.data_file;

                const Field<T>* fp = nullptr;
                if (!token.empty()) {
                    FieldToken ft = parse_field_token(token);
                    auto mat = get_matrix(ft.path);           // shared, loaded once
                    extract_column(*mat, ft.col, local);
                    fp = &local;
                }

                ViewResult<T> r = process_view<T>({(T)b.nx, (T)b.ny, (T)b.nz}, mesh, fp,
                                                  view_res[k], cull, renderer);
                if (!out_pre.empty()) {
                    std::ostringstream oss; oss << out_pre << "_" << std::setw(4) << std::setfill('0') << k << ".ppm";
                    r.output_image = oss.str();
                    save_ppm(renderer, r.output_image, r.min_val, r.max_val, r.has_field);
                }
                results[k] = std::move(r);
            } catch (const std::exception& e) {
                std::lock_guard<std::mutex> lk(err_mtx);
                if (err_msg.empty()) err_msg = "timestep " + std::to_string(k) + ": " + e.what();
                next.store(batch.size());   // signal other workers to stop
                return;
            }
        }
    };

    unsigned n_workers = std::max(1u, std::min<unsigned>(threads, (unsigned)batch.size()));
    if (n_workers <= 1) {
        worker();
    } else {
        std::vector<std::thread> pool;
        for (unsigned t=0; t<n_workers; ++t) pool.emplace_back(worker);
        for (auto& th : pool) th.join();
    }
    if (!err_msg.empty()) throw std::runtime_error(err_msg);

    // Emit results in view order (after compute, so threaded output stays ordered).
    if (json) {
        std::cout << "{\n  \"meta\": {\n"
                  << "    \"mesh\": \"" << esc(obj) << "\",\n"
                  << "    \"precision\": \"" << (sizeof(T)==4?"float":"double") << "\",\n"
                  << "    \"views\": " << batch.size() << "\n"
                  << "  },\n  \"results\": [\n";
    }
    for (size_t k=0; k<batch.size(); ++k) {
        const auto& b = batch[k];
        const auto& r = results[k];
        if (json) {
            std::cout << "    {\n"
                      << "      \"idx\": " << k << ",\n"
                      << "      \"normal\": [" << jnum(b.nx) << "," << jnum(b.ny) << "," << jnum(b.nz) << "],\n"
                      << "      \"resolution\": " << jnum(view_res[k]) << ",\n"
                      << "      \"area\": " << jnum(r.area) << ",\n"
                      << "      \"pixels\": " << r.covered_pixels << ",\n"
                      << "      \"width\": " << r.image_width << ",\n"
                      << "      \"height\": " << r.image_height << ",\n";
            if (r.has_field) {
                std::cout << "      \"average\": " << jnum(r.average_value) << ",\n"
                          << "      \"integral\": " << jnum(r.integral) << ",\n"
                          << "      \"min\": " << jnum(r.min_val) << ",\n"
                          << "      \"max\": " << jnum(r.max_val) << ",\n";
            } else {
                std::cout << "      \"average\": null,\n"
                          << "      \"integral\": null,\n"
                          << "      \"min\": null,\n"
                          << "      \"max\": null,\n";
            }
            std::cout << "      \"image\": " << (r.output_image.empty() ? "null" : "\"" + esc(r.output_image) + "\"") << ",\n"
                      << "      \"time\": " << jnum(r.time_seconds) << "\n"
                      << "    }" << (k==batch.size()-1?"":",") << "\n";
        } else {
            std::cout << "View " << k << " | Area: " << r.area;
            if (r.has_field)
                std::cout << " | Avg: " << r.average_value << " | Integral: " << r.integral;
            std::cout << " | Pixels: " << r.covered_pixels
                      << " | Time: " << r.time_seconds << "s\n";
        }
    }
    if (json) std::cout << "  ]\n}\n";
}

template void run_app<float>(const std::string&, const std::string&, const std::string&,
                             const std::vector<BatchEntry>&, double, bool, bool, unsigned, ValueMode);
template void run_app<double>(const std::string&, const std::string&, const std::string&,
                              const std::vector<BatchEntry>&, double, bool, bool, unsigned, ValueMode);

} // namespace flatland
