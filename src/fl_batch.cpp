#include "fl_batch.hpp"

#include "fl_field.hpp"
#include "fl_image.hpp"
#include "fl_io.hpp"
#include "fl_parallel.hpp"
#include "fl_project.hpp"
#include "fl_raster.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <thread>
#include <vector>
#include <cstddef>
#include <cstdio>

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

namespace {

// Rough in-memory size of a loaded matrix, for the cache budget below.
template <typename T>
size_t matrix_bytes(const FieldMatrix<T>& m) { return m.data.size() * sizeof(T) + sizeof(m); }

// Keep at most this much field data resident. A batch that names one file per
// timestep would otherwise grow the cache without limit — roughly 16 MB per
// distinct file on a 10k-vertex mesh with 400 columns, which a long series will
// turn into an out-of-memory kill.
const size_t CACHE_BUDGET_BYTES = 256u << 20;   // 256 MB

} // namespace

/* ----------------------
   Application Runner
   ---------------------- */
template <typename T>
void run_app(const std::string& obj, const std::string& default_data, const std::string& out_pre,
             const std::vector<BatchEntry>& batch, double default_res, bool cull, bool json,
             unsigned threads, ValueMode forced_mode) {

    Mesh<T> mesh = load_mesh<T>(obj);

    // Field-matrix cache.
    //
    // Each distinct file is parsed exactly once and shared read-only; workers then
    // extract their own timestep's column from it concurrently. Loading happens
    // OUTSIDE the map lock, via a per-entry std::call_once, so two workers wanting
    // two different files do not serialize behind each other — only two workers
    // wanting the SAME file wait, which is the point. A load that fails is
    // remembered, so a bad path is not re-parsed by every worker in turn.
    //
    // The cache is bounded. When it exceeds its budget, entries no one is
    // currently holding are dropped; anything in use is kept, so eviction can
    // never pull a matrix out from under a worker.
    struct CacheSlot {
        std::once_flag once;
        std::shared_ptr<FieldMatrix<T>> mat;
        std::string error;
        size_t bytes = 0;
    };
    std::map<std::string, std::shared_ptr<CacheSlot>> mat_cache;
    std::mutex cache_mtx;

    auto get_matrix = [&](const std::string& path) -> std::shared_ptr<FieldMatrix<T>> {
        std::shared_ptr<CacheSlot> slot;
        {
            std::lock_guard<std::mutex> lk(cache_mtx);
            auto& s = mat_cache[path];
            if (!s) s = std::make_shared<CacheSlot>();
            slot = s;
        }
        std::call_once(slot->once, [&]() {
            try {
                auto m = std::make_shared<FieldMatrix<T>>();
                load_matrix_into(path, mesh, *m, forced_mode);
                slot->mat = std::move(m);
            } catch (const std::exception& e) {
                slot->error = e.what();     // memoized: do not re-parse a bad file
            }
        });
        if (!slot->mat)
            throw std::runtime_error(slot->error.empty()
                                     ? "cannot load data file '" + path + "'" : slot->error);

        {
            std::lock_guard<std::mutex> lk(cache_mtx);
            // `bytes` is accounted here rather than inside call_once: the
            // eviction scan below reads every slot's size while holding this
            // lock, so writing it unlocked from a loading thread is a race.
            // call_once already established the happens-before for slot->mat.
            if (slot->bytes == 0) slot->bytes = matrix_bytes(*slot->mat);
            size_t total = 0;
            for (const auto& kv : mat_cache) total += kv.second->bytes;
            for (auto it = mat_cache.begin(); it != mat_cache.end() && total > CACHE_BUDGET_BYTES; ) {
                // use_count()==1 means only the map holds this slot, so no worker
                // is reading it and dropping it is safe. A slot evicted while its
                // file is still needed simply reloads on the next request.
                if (it->second != slot && it->second.use_count() == 1) {
                    total -= std::min(total, it->second->bytes);
                    it = mat_cache.erase(it);
                } else {
                    ++it;
                }
            }
        }
        return slot->mat;
    };

    // Pre-warm + validate the default field up front so a bad -d fails fast.
    if (!default_data.empty()) get_matrix(parse_field_token(default_data).path);

    std::vector<T> view_res(batch.size());
    for (size_t k=0; k<batch.size(); ++k)
        view_res[k] = (batch[k].resolution > 0) ? (T)batch[k].resolution : (T)default_res;

    std::vector<ViewResult<T>> results(batch.size());

    // Images are written as they are produced, so peak memory stays O(threads)
    // rather than O(timesteps). If the batch then fails we remove what was
    // written, rather than leaving a partial series on disk for a run that
    // reported nothing.
    std::mutex img_mtx;
    std::vector<std::string> written_images;

    const unsigned n_workers = choose_workers(batch.size(), threads);
    if (!json)
        std::cerr << "Running " << batch.size() << " view(s) with "
                  << (sizeof(T)==4 ? "float" : "double") << " precision on "
                  << n_workers << " thread(s)...\n";

    // Per-worker scratch, indexed by the id parallel_for hands each worker.
    // Keeping it here rather than in thread_local storage means it is released
    // when this call returns, and makes the ownership obvious.
    std::vector<Renderer<T>> renderers(n_workers);
    std::vector<Field<T>> locals(n_workers);

    try {
        parallel_for(batch.size(), n_workers, [&](size_t k, unsigned w) {
            const auto& b = batch[k];
            const std::string& token = b.data_file.empty() ? default_data : b.data_file;

            const Field<T>* fp = nullptr;
            if (!token.empty()) {
                FieldToken ft = parse_field_token(token);
                auto mat = get_matrix(ft.path);           // shared, loaded once
                extract_column(*mat, ft.col, locals[w]);
                fp = &locals[w];
            }

            ViewResult<T> r = process_view<T>({(T)b.nx, (T)b.ny, (T)b.nz}, mesh, fp,
                                              view_res[k], cull, renderers[w]);

            // A view that covered nothing has no raster to write. Writing one
            // anyway would emit the previous view's image, or a malformed 0x0
            // file; the JSON reports image: null instead.
            if (!out_pre.empty() && r.image_width > 0 && r.image_height > 0) {
                std::ostringstream oss;
                oss << out_pre << "_" << std::setw(4) << std::setfill('0') << k << ".ppm";
                r.output_image = oss.str();
                save_ppm(renderers[w], r.output_image, r.min_val, r.max_val, r.has_stats);
                std::lock_guard<std::mutex> lk(img_mtx);
                written_images.push_back(r.output_image);
            }
            results[k] = std::move(r);
        });
    } catch (const ParallelError& e) {
        for (const auto& p : written_images) std::remove(p.c_str());
        throw std::runtime_error("timestep " + std::to_string(e.index) + ": " + e.what());
    }

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
            // has_stats, not has_field: a view carrying a field but covering no
            // pixels has no statistics, and emitting 0 would be a fabrication a
            // consumer could not distinguish from a measurement.
            if (r.has_stats) {
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
            if (r.has_stats)
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
