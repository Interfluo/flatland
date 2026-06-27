//
//  main.cpp
//  FlatLand
//
//  Orthographic projection + rasterization for geometric analysis of 3D meshes.
//
//  Core capabilities:
//   - Projected (silhouette / visible) surface area from arbitrary view directions
//   - Generic per-view field statistics over the visible projection:
//       average, min, max, and the AREA INTEGRAL  ∫ f dA = Σ value · pixel_area.
//     The integral is the building block for quantities such as radiant
//     intensity (supply a radiance field) — FlatLand stays domain-agnostic.
//   - Runtime precision switching (float / double)
//   - Batch mode with per-view resolution and data files
//   - Parallel batch processing (standard library threads only)
//   - PPM heatmap export and structured JSON output
//
//  No third-party dependencies. C++17.
//

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <memory>
#include <thread>
#include <atomic>
#include <mutex>
#include <cstdlib>
#include <clocale>
#include <stdexcept>

/* ----------------------
   Math & Structs (Templated)
   ---------------------- */
template <typename T> struct Vec3 { T x, y, z; };
template <typename T> struct Vec2 { T x, y; };

template <typename T>
struct TriProjected {
    Vec2<T> v0, v1, v2;
    T z0, z1, z2;
    T val0, val1, val2;
    T inv_area;
};

struct Triangle3D {
    int v0_idx, v1_idx, v2_idx; // Indices into the vertex array
};

enum ValueMode { MODE_NONE, MODE_NODE, MODE_FACE };

// Read-only geometry, loaded once and shared across all views/threads.
template <typename T>
struct Mesh {
    std::vector<Vec3<T>> vertices;
    std::vector<Triangle3D> faces;
};

// A scalar field defined on the mesh. Kept separate from Mesh so different
// views can use different fields without mutating shared state (thread-safe).
template <typename T>
struct Field {
    std::vector<T> values;
    ValueMode mode = MODE_NONE;
};

template <typename T>
struct ViewResult {
    bool has_field = false;
    T area = 0;
    T average_value = 0;
    T integral = 0;     // Σ value · pixel_area  (the area integral of the field)
    T min_val = 0;
    T max_val = 0;
    double time_seconds = 0.0;
    long covered_pixels = 0;
    int image_width = 0;
    int image_height = 0;
    std::string output_image;
};

struct BatchEntry {
    double nx, ny, nz;
    double resolution = -1.0; // -1 means use global default
    std::string data_file;    // empty means use global default
};

/* ----------------------
   Math Ops
   ---------------------- */
template <typename T> inline T dot(const Vec3<T>& a, const Vec3<T>& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
template <typename T> inline Vec3<T> cross(const Vec3<T>& a, const Vec3<T>& b) { return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x}; }
template <typename T> inline Vec3<T> operator-(const Vec3<T>& a, const Vec3<T>& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
template <typename T> inline Vec3<T> normalize(const Vec3<T>& v) { T n = std::sqrt(dot(v,v)); return n>0 ? Vec3<T>{v.x/n, v.y/n, v.z/n} : Vec3<T>{0,0,0}; }

/* ----------------------
   Renderer (one per worker thread)
   ---------------------- */
template <typename T>
class Renderer {
public:
    int Nx = 0, Ny = 0;
    std::vector<T> zbuffer;
    std::vector<T> val_buffer;
    std::vector<uint8_t> mask;

    void resize(int w, int h, bool use_vals) {
        size_t needed = (size_t)w * h;
        if (needed == 0) { Nx = Ny = 0; return; }
        Nx = w; Ny = h;
        if (zbuffer.size() < needed) zbuffer.resize(needed);
        if (mask.size() < needed) mask.resize(needed);
        if (use_vals && val_buffer.size() < needed) val_buffer.resize(needed);

        std::fill(zbuffer.begin(), zbuffer.begin() + needed, std::numeric_limits<T>::max());
        std::fill(mask.begin(), mask.begin() + needed, 0);
        if (use_vals) std::fill(val_buffer.begin(), val_buffer.begin() + needed, 0);
    }
};

/* ----------------------
   Rasterization Logic
   ---------------------- */
template <typename T>
inline T edge_eval(const Vec2<T>& a, const Vec2<T>& b, const Vec2<T>& p) {
    return (p.x - a.x)*(b.y - a.y) - (p.y - a.y)*(b.x - a.x);
}

template <typename T>
void rasterize(const TriProjected<T>& tri, Renderer<T>& r, T xmin, T ymin, T pix_sz, bool use_vals) {
    T min_x = std::min({tri.v0.x, tri.v1.x, tri.v2.x});
    T max_x = std::max({tri.v0.x, tri.v1.x, tri.v2.x});
    T min_y = std::min({tri.v0.y, tri.v1.y, tri.v2.y});
    T max_y = std::max({tri.v0.y, tri.v1.y, tri.v2.y});

    int ix_min = std::max(0, static_cast<int>((min_x - xmin)/pix_sz));
    int ix_max = std::min(r.Nx - 1, static_cast<int>((max_x - xmin)/pix_sz));
    int iy_min = std::max(0, static_cast<int>((min_y - ymin)/pix_sz));
    int iy_max = std::min(r.Ny - 1, static_cast<int>((max_y - ymin)/pix_sz));

    Vec2<T> p_start = { xmin + (ix_min + (T)0.5)*pix_sz, ymin + (iy_min + (T)0.5)*pix_sz };

    // Edge function values at the first sampled pixel center
    T w0_row = edge_eval(tri.v1, tri.v2, p_start);
    T w1_row = edge_eval(tri.v2, tri.v0, p_start);
    T w2_row = edge_eval(tri.v0, tri.v1, p_start);

    // Per-pixel increments
    T A0 = (tri.v2.y - tri.v1.y)*pix_sz, B0 = (tri.v1.x - tri.v2.x)*pix_sz;
    T A1 = (tri.v0.y - tri.v2.y)*pix_sz, B1 = (tri.v2.x - tri.v0.x)*pix_sz;
    T A2 = (tri.v1.y - tri.v0.y)*pix_sz, B2 = (tri.v0.x - tri.v1.x)*pix_sz;

    for (int iy = iy_min; iy <= iy_max; ++iy) {
        T w0=w0_row, w1=w1_row, w2=w2_row;
        size_t idx = (size_t)iy * r.Nx + ix_min;
        for (int ix = ix_min; ix <= ix_max; ++ix) {
            if (w0 >= 0 && w1 >= 0 && w2 >= 0) {
                T z = (w0*tri.inv_area)*tri.z0 + (w1*tri.inv_area)*tri.z1 + (w2*tri.inv_area)*tri.z2;
                if (z < r.zbuffer[idx]) {
                    r.zbuffer[idx] = z;
                    r.mask[idx] = 1;
                    if (use_vals)
                        r.val_buffer[idx] = (w0*tri.val0 + w1*tri.val1 + w2*tri.val2)*tri.inv_area;
                }
            }
            w0+=A0; w1+=A1; w2+=A2; idx++;
        }
        w0_row+=B0; w1_row+=B1; w2_row+=B2;
    }
}

/* ----------------------
   IO & Loaders
   ---------------------- */
// Templated OBJ loader to support float/double geometry.
template <typename T>
Mesh<T> load_obj(const std::string& filename) {
    Mesh<T> mesh;
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("cannot open mesh file '" + filename + "'");

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (line[0] == 'v' && line[1] == ' ') {
            double x,y,z;
            if(sscanf(line.c_str(), "v %lf %lf %lf", &x, &y, &z) == 3)
                mesh.vertices.push_back({(T)x, (T)y, (T)z});
        } else if (line[0] == 'f') {
            std::vector<int> idxs;
            const char* s = line.c_str() + 1;
            while (*s) {
                while (*s == ' ' || *s == '\t') s++;
                if (!*s) break;
                int idx = std::atoi(s);
                idxs.push_back(idx > 0 ? idx - 1 : idx + (int)mesh.vertices.size());
                while (*s && *s != ' ' && *s != '\t') s++;
            }
            for (size_t k=2; k<idxs.size(); ++k)
                mesh.faces.push_back({idxs[0], idxs[k-1], idxs[k]});
        }
    }
    if (mesh.vertices.empty()) throw std::runtime_error("mesh '" + filename + "' contains no vertices");
    if (mesh.faces.empty())    throw std::runtime_error("mesh '" + filename + "' contains no faces");

    // Guard against malformed files: every face index must reference a real vertex,
    // otherwise process_view would read out of bounds.
    const int nv = (int)mesh.vertices.size();
    for (const auto& f : mesh.faces) {
        if (f.v0_idx < 0 || f.v0_idx >= nv ||
            f.v1_idx < 0 || f.v1_idx >= nv ||
            f.v2_idx < 0 || f.v2_idx >= nv)
            throw std::runtime_error("mesh '" + filename + "' references an out-of-range vertex index");
    }
    return mesh;
}

// Fast field loader: slurp the whole file and parse with strtod (far faster than
// iostream extraction for the time-series case of thousands of field files).
// Reuses `field.values` storage so per-timestep reloads don't churn allocations.
template <typename T>
void load_field_into(const std::string& filename, const Mesh<T>& mesh, Field<T>& field) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("cannot open data file '" + filename + "'");
    std::string buf((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

    field.values.clear();
    const char* p = buf.c_str();
    char* end;
    double d = std::strtod(p, &end);
    while (end != p) {                 // strtod skips leading whitespace itself
        field.values.push_back((T)d);
        p = end;
        d = std::strtod(p, &end);
    }

    if (field.values.size() == mesh.vertices.size())   field.mode = MODE_NODE;
    else if (field.values.size() == mesh.faces.size()) field.mode = MODE_FACE;
    else throw std::runtime_error("data file '" + filename + "' has " +
            std::to_string(field.values.size()) + " values, but mesh has " +
            std::to_string(mesh.vertices.size()) + " vertices / " +
            std::to_string(mesh.faces.size()) + " faces");
}

/* ----------------------
   Image Output
   ---------------------- */
template <typename T>
void save_ppm(const Renderer<T>& r, const std::string& fname, T min_v, T max_v, bool use_vals) {
    std::ofstream f(fname, std::ios::binary);
    if (!f.is_open()) throw std::runtime_error("cannot write image '" + fname + "'");
    f << "P6\n" << r.Nx << " " << r.Ny << "\n255\n";
    double range = (double)(max_v - min_v);
    if (std::abs(range) < 1e-9) range = 1.0;
    const uint8_t bg[3] = {30,30,35};

    // Simple 5-stop heatmap (blue -> cyan -> yellow -> orange -> red)
    const double stops[5][3] = {{0,0,1}, {0,1,1}, {1,1,0}, {1,.5,0}, {1,0,0}};

    for (int iy=r.Ny-1; iy>=0; --iy) {
        for (int ix=0; ix<r.Nx; ++ix) {
            size_t i = (size_t)iy * r.Nx + ix;
            if (!r.mask[i]) { f.write((const char*)bg, 3); continue; }
            uint8_t rgb[3];
            if (!use_vals) {
                rgb[0]=rgb[1]=rgb[2]=255;            // silhouette: solid white
            } else {
                double t = ((double)r.val_buffer[i] - min_v)/range;
                t = std::max(0.0, std::min(1.0, t)) * 4.0;
                int idx = (int)t;
                if (idx >= 4) idx = 3;
                double fr = t - idx;
                for(int c=0;c<3;++c) rgb[c] = (uint8_t)(((1.0-fr)*stops[idx][c] + fr*stops[idx+1][c])*255);
            }
            f.write((const char*)rgb, 3);
        }
    }
}

/* ----------------------
   Processing Pipeline
   ---------------------- */
template <typename T>
ViewResult<T> process_view(const Vec3<T>& view_dir, const Mesh<T>& mesh, const Field<T>* field,
                           T pix_sz, bool cull, Renderer<T>& renderer) {
    auto t0 = std::chrono::high_resolution_clock::now();
    ViewResult<T> res;
    const bool use_vals = (field && field->mode != MODE_NONE);
    res.has_field = use_vals;

    // Orthonormal view basis (n = view normal, u/v = screen axes)
    Vec3<T> n = normalize(view_dir);
    Vec3<T> h = (std::abs(n.x) < 0.9) ? Vec3<T>{1,0,0} : Vec3<T>{0,1,0};
    Vec3<T> u = normalize(cross(h, n));
    Vec3<T> v = cross(n, u);

    // Project triangles to the view plane
    std::vector<TriProjected<T>> tris;
    tris.reserve(mesh.faces.size());
    T b_min_x=1e20, b_max_x=-1e20, b_min_y=1e20, b_max_y=-1e20;

    for (size_t i=0; i<mesh.faces.size(); ++i) {
        const auto& f = mesh.faces[i];
        const Vec3<T>& v0 = mesh.vertices[f.v0_idx];
        const Vec3<T>& v1 = mesh.vertices[f.v1_idx];
        const Vec3<T>& v2 = mesh.vertices[f.v2_idx];

        Vec3<T> tri_n = cross(v1-v0, v2-v0);
        if (cull && dot(tri_n, n) <= 0) continue;   // keep faces pointing toward the viewer

        Vec2<T> p0 = {dot(v0,u), dot(v0,v)};
        Vec2<T> p1 = {dot(v1,u), dot(v1,v)};
        Vec2<T> p2 = {dot(v2,u), dot(v2,v)};

        T val0=0, val1=0, val2=0;
        if (use_vals) {
            if (field->mode == MODE_NODE) {
                val0=field->values[f.v0_idx]; val1=field->values[f.v1_idx]; val2=field->values[f.v2_idx];
            } else { // MODE_FACE
                val0=val1=val2=field->values[i];
            }
        }

        T area2 = edge_eval(p0, p1, p2); // twice signed screen-space area; >0 = CCW
        if (std::abs(area2) < 1e-12) continue;

        b_min_x = std::min({b_min_x, p0.x, p1.x, p2.x});
        b_max_x = std::max({b_max_x, p0.x, p1.x, p2.x});
        b_min_y = std::min({b_min_y, p0.y, p1.y, p2.y});
        b_max_y = std::max({b_max_y, p0.y, p1.y, p2.y});

        // Normalize winding to CCW so edge functions are positive inside.
        if (area2 < 0) {
            tris.push_back({ p0, p2, p1,
                             dot(v0,n), dot(v2,n), dot(v1,n),
                             val0, val2, val1,
                             (T)(1.0/(-area2)) });
        } else {
            tris.push_back({ p0, p1, p2,
                             dot(v0,n), dot(v1,n), dot(v2,n),
                             val0, val1, val2,
                             (T)(1.0/area2) });
        }
    }

    if (!tris.empty()) {
        b_min_x -= pix_sz; b_max_x += pix_sz;
        b_min_y -= pix_sz; b_max_y += pix_sz;
        int w = (int)std::ceil((b_max_x - b_min_x)/pix_sz);
        int h = (int)std::ceil((b_max_y - b_min_y)/pix_sz);

        renderer.resize(w, h, use_vals);
        res.image_width=w; res.image_height=h;

        for (const auto& t : tris) rasterize(t, renderer, b_min_x, b_min_y, pix_sz, use_vals);

        // Aggregate statistics over visible pixels
        size_t npix = (size_t)w * h;
        double sum = 0;
        T mn = std::numeric_limits<T>::max();
        T mx = -std::numeric_limits<T>::max();
        for (size_t k=0; k<npix; ++k) {
            if (!renderer.mask[k]) continue;
            res.covered_pixels++;
            if (use_vals) {
                T val = renderer.val_buffer[k];
                sum += val;
                if (val < mn) mn = val;
                if (val > mx) mx = val;
            }
        }
        T pix_area = pix_sz * pix_sz;
        res.area = res.covered_pixels * pix_area;
        if (use_vals && res.covered_pixels) {
            res.average_value = (T)(sum / res.covered_pixels);
            res.integral = (T)(sum * pix_area);  // ∫ f dA over the visible projection
            res.min_val = mn;
            res.max_val = mx;
        }
    }

    res.time_seconds = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0).count();
    return res;
}

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
             unsigned threads) {

    Mesh<T> mesh = load_obj<T>(obj);

    // A fixed field (-d) is loaded ONCE and shared read-only across all timesteps/threads.
    std::shared_ptr<Field<T>> default_field;
    if (!default_data.empty()) {
        default_field = std::make_shared<Field<T>>();
        load_field_into(default_data, mesh, *default_field);
    }

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
        Field<T> local;            // scratch for per-view data files
        std::string local_path;    // last file loaded into `local`
        size_t k;
        while ((k = next.fetch_add(1)) < batch.size()) {
            try {
                const auto& b = batch[k];
                const std::string& path = b.data_file.empty() ? default_data : b.data_file;

                const Field<T>* fp = nullptr;
                if (path.empty()) {
                    fp = nullptr;
                } else if (path == default_data) {
                    fp = default_field.get();                 // shared, already loaded
                } else {
                    if (path != local_path) { load_field_into(path, mesh, local); local_path = path; }
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

/* ----------------------
   CLI Parsing
   ---------------------- */
void print_help(const char* prog) {
    std::cout << "FlatLand: Mesh Projection & Rasterization Tool\n"
              << "Usage: " << prog << " <mesh.obj> [options]\n\n"
              << "Core Options:\n"
              << "  -v, --view <x> <y> <z>   Add a single view direction. Repeatable.\n"
              << "  -b, --batch <file>       Load views from a file (see 'Batch Format').\n"
              << "                           Views from -v and -b are additive.\n\n"
              << "Analysis Parameters:\n"
              << "  -r, --res <val>          Default pixel resolution (default: 0.001).\n"
              << "                           Used when a batch entry omits its own resolution.\n"
              << "  -d, --data <file>        Default scalar field file (one value per line;\n"
              << "                           length must equal the vertex or face count).\n"
              << "  -p, --precision <mode>   Math precision: 'float' (default) or 'double'.\n"
              << "      --no-cull            Disable backface culling (render all faces).\n"
              << "  -t, --threads <n>        Worker threads for batch views (default: auto).\n\n"
              << "Output Control:\n"
              << "  -o, --out <prefix>       Save PPM heatmaps as <prefix>_<idx>.ppm.\n"
              << "  -j, --json               Emit structured JSON to stdout.\n"
              << "  -h, --help               Show this help message.\n\n"
              << "Per-View Outputs:\n"
              << "  area      visible projected area (covered_pixels * resolution^2)\n"
              << "  average   mean field value over visible pixels (with -d / batch data)\n"
              << "  integral  area integral of the field: sum(value * pixel_area).\n"
              << "            Supply a radiance field to get radiant intensity, etc.\n"
              << "  min/max   field extremes over visible pixels\n\n"
              << "Batch File Format:\n"
              << "  One view per line; '#' starts a comment.\n"
              << "  Columns: <nx> <ny> <nz> [resolution] [data_file_path]\n"
              << "    1.0 0.0 0.0                 (default resolution & data)\n"
              << "    0.0 1.0 0.0 0.005           (override resolution)\n"
              << "    0.0 0.0 1.0 0.002 data.txt  (override both)\n";
}

// Parse the value following a flag; errors if it is missing.
static std::string need_value(int& i, int argc, char* argv[], const std::string& flag) {
    if (i+1 >= argc) throw std::runtime_error("option '" + flag + "' requires an argument");
    return argv[++i];
}

static double parse_double(const std::string& s, const std::string& ctx) {
    try {
        size_t pos; double d = std::stod(s, &pos);
        if (pos != s.size()) throw std::invalid_argument("");
        return d;
    } catch (...) { throw std::runtime_error("invalid number '" + s + "' for " + ctx); }
}

int main(int argc, char* argv[]) {
    // Force the C locale so numeric parsing (strtod/sscanf) always uses '.' as the
    // decimal separator, regardless of the user's environment locale.
    std::setlocale(LC_ALL, "C");

    if (argc < 2) { print_help(argv[0]); return 1; }

    std::string obj_file, data_file, batch_file, out_pre, prec="float";
    double res = 0.001;
    bool json=false, cull=true;
    unsigned threads = 0; // 0 => auto
    std::vector<BatchEntry> batch;

    try {
        for (int i=1; i<argc; ++i) {
            std::string a = argv[i];
            if      (a == "-h" || a == "--help")      { print_help(argv[0]); return 0; }
            else if (a == "-d" || a == "--data")      { data_file  = need_value(i,argc,argv,a); }
            else if (a == "-b" || a == "--batch")     { batch_file = need_value(i,argc,argv,a); }
            else if (a == "-o" || a == "--out")       { out_pre    = need_value(i,argc,argv,a); }
            else if (a == "-r" || a == "--res")       { res = parse_double(need_value(i,argc,argv,a), "-r/--res"); }
            else if (a == "-p" || a == "--precision") { prec = need_value(i,argc,argv,a); }
            else if (a == "-t" || a == "--threads")   { threads = (unsigned)parse_double(need_value(i,argc,argv,a), "-t/--threads"); }
            else if (a == "-j" || a == "--json")      { json = true; }
            else if (a == "--no-cull")                { cull = false; }
            else if (a == "-v" || a == "--view") {
                if (i+3 >= argc) throw std::runtime_error("option '" + a + "' requires three numbers: x y z");
                double x = parse_double(argv[i+1], "-v x");
                double y = parse_double(argv[i+2], "-v y");
                double z = parse_double(argv[i+3], "-v z");
                if (x==0 && y==0 && z==0) throw std::runtime_error("view direction cannot be the zero vector");
                batch.push_back({x, y, z, -1.0, ""});
                i += 3;
            }
            else if (!a.empty() && a[0] == '-') {
                throw std::runtime_error("unknown option '" + a + "' (try --help)");
            }
            else if (obj_file.empty()) { obj_file = a; }
            else throw std::runtime_error("unexpected argument '" + a + "' (mesh already set to '" + obj_file + "')");
        }

        if (obj_file.empty()) throw std::runtime_error("no mesh (.obj) file specified");
        if (prec != "float" && prec != "double") throw std::runtime_error("precision must be 'float' or 'double', got '" + prec + "'");
        if (res <= 0) throw std::runtime_error("resolution must be positive");

        // Parse batch file (additive with any -v views)
        if (!batch_file.empty()) {
            std::ifstream f(batch_file);
            if (!f.is_open()) throw std::runtime_error("cannot open batch file '" + batch_file + "'");
            // Relative data-file paths in a batch are resolved against the batch file's
            // directory, so a committed time-series case works from any CWD.
            std::string batch_dir;
            size_t slash = batch_file.find_last_of('/');
            if (slash != std::string::npos) batch_dir = batch_file.substr(0, slash + 1);
            std::string line;
            int lineno = 0;
            while (std::getline(f, line)) {
                ++lineno;
                size_t c = line.find('#');
                if (c != std::string::npos) line = line.substr(0, c);
                std::stringstream ss(line);
                double nx, ny, nz;
                if (ss >> nx >> ny >> nz) {
                    if (nx==0 && ny==0 && nz==0)
                        throw std::runtime_error("batch line " + std::to_string(lineno) + ": zero view vector");
                    BatchEntry be = {nx, ny, nz, -1.0, ""};
                    auto resolve = [&](const std::string& p) {
                        return (!p.empty() && p[0] == '/') ? p : batch_dir + p;
                    };
                    // The optional trailing tokens are [resolution] and/or [data_file],
                    // in that order. We sniff each token: a pure number is a resolution,
                    // anything else is a data-file path. This lets a view attach data
                    // without having to restate the default resolution.
                    std::string tok;
                    if (ss >> tok) {
                        char* endp; double rv = std::strtod(tok.c_str(), &endp);
                        if (*endp == '\0') { // pure number -> resolution
                            if (rv <= 0) throw std::runtime_error("batch line " + std::to_string(lineno) + ": resolution must be positive");
                            be.resolution = rv;
                            std::string d_temp;
                            if (ss >> d_temp) be.data_file = resolve(d_temp);
                        } else {             // non-numeric -> data file, default resolution
                            be.data_file = resolve(tok);
                        }
                    }
                    batch.push_back(be);
                } else if (!line.empty() && line.find_first_not_of(" \t\r\n") != std::string::npos) {
                    throw std::runtime_error("batch line " + std::to_string(lineno) + ": expected '<nx> <ny> <nz> ...'");
                }
            }
        }

        if (batch.empty()) throw std::runtime_error("no views specified (use -v or -b)");

        unsigned hw = std::thread::hardware_concurrency();
        if (threads == 0) threads = hw ? hw : 1;

        if (!json) std::cerr << "Running " << batch.size() << " view(s) with " << prec
                             << " precision on " << threads << " thread(s)...\n";

        if (prec == "float") run_app<float>(obj_file, data_file, out_pre, batch, res, cull, json, threads);
        else                 run_app<double>(obj_file, data_file, out_pre, batch, res, cull, json, threads);

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
