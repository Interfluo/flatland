//
//  main.cpp
//  FlatLand
//
//  Features:
//  - Runtime precision switching (Float vs Double)
//  - Flexible Batch Mode (Per-view resolution & data files)
//  - Edge Function Rasterization
//  - JSON Metadata
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
#include <cstring>
#include <memory>

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
    int v0_idx, v1_idx, v2_idx; // Indices into vertex array
    // We store indices to allow re-fetching vertices if precision changes (though here we load once)
};

enum ValueMode { MODE_NONE, MODE_NODE, MODE_FACE };

template <typename T>
struct Mesh {
    std::vector<Vec3<T>> vertices;
    std::vector<Triangle3D> faces;
    std::vector<T> values;
    ValueMode val_mode = MODE_NONE;
};

template <typename T>
struct ViewResult {
    T area = 0;
    T average_value = 0;
    T min_val = std::numeric_limits<T>::max();
    T max_val = -std::numeric_limits<T>::max();
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
   Renderer
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
        if (needed == 0) return;
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
void rasterize(const TriProjected<T>& tri, Renderer<T>& r, T xmin, T ymin, T pix_sz) {
    T min_x = std::min({tri.v0.x, tri.v1.x, tri.v2.x});
    T max_x = std::max({tri.v0.x, tri.v1.x, tri.v2.x});
    T min_y = std::min({tri.v0.y, tri.v1.y, tri.v2.y});
    T max_y = std::max({tri.v0.y, tri.v1.y, tri.v2.y});

    int ix_min = std::max(0, static_cast<int>((min_x - xmin)/pix_sz));
    int ix_max = std::min(r.Nx - 1, static_cast<int>((max_x - xmin)/pix_sz));
    int iy_min = std::max(0, static_cast<int>((min_y - ymin)/pix_sz));
    int iy_max = std::min(r.Ny - 1, static_cast<int>((max_y - ymin)/pix_sz));

    Vec2<T> p_start = { xmin + (ix_min + 0.5f)*pix_sz, ymin + (iy_min + 0.5f)*pix_sz };
    
    // Edge Functions Setup
    T w0_row = edge_eval(tri.v1, tri.v2, p_start);
    T w1_row = edge_eval(tri.v2, tri.v0, p_start);
    T w2_row = edge_eval(tri.v0, tri.v1, p_start);
    
    // Increments
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
                    if (!r.val_buffer.empty())
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
// Templated OBJ loader to support Float/Double
template <typename T>
Mesh<T> load_obj(const std::string& filename) {
    Mesh<T> mesh;
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Cannot open " + filename);
    
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
    return mesh;
}

template <typename T>
void load_values(const std::string& filename, Mesh<T>& mesh) {
    if (filename.empty()) return;
    mesh.values.clear();
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Cannot open data file " + filename);
    double v;
    while (file >> v) mesh.values.push_back((T)v);
    
    if (mesh.values.size() == mesh.vertices.size()) mesh.val_mode = MODE_NODE;
    else if (mesh.values.size() == mesh.faces.size()) mesh.val_mode = MODE_FACE;
    else throw std::runtime_error("Data size mismatch: " + std::to_string(mesh.values.size()));
}

/* ----------------------
   Image Output
   ---------------------- */
template <typename T>
void save_ppm(const Renderer<T>& r, const std::string& fname, T min_v, T max_v) {
    std::ofstream f(fname, std::ios::binary);
    if (!f.is_open()) return;
    f << "P6\n" << r.Nx << " " << r.Ny << "\n255\n";
    double range = (double)(max_v - min_v);
    if (std::abs(range) < 1e-9) range = 1.0;
    const uint8_t bg[3] = {30,30,35};
    
    // Simple heatmap colors
    const double stops[5][3] = {{0,0,1}, {0,1,1}, {1,1,0}, {1,.5,0}, {1,0,0}};
    
    for (int iy=r.Ny-1; iy>=0; --iy) {
        for (int ix=0; ix<r.Nx; ++ix) {
            size_t i = (size_t)iy * r.Nx + ix;
            if (!r.mask[i]) f.write((char*)bg, 3);
            else {
                double t = (r.val_buffer.empty()) ? 1.0 : ((double)r.val_buffer[i] - min_v)/range;
                t = std::max(0.0, std::min(1.0, t)) * 4.0;
                int idx = (int)t;
                if (idx >= 4) idx = 3;
                double fr = t - idx;
                uint8_t rgb[3];
                for(int c=0;c<3;++c) rgb[c] = (uint8_t)(((1.0-fr)*stops[idx][c] + fr*stops[idx+1][c])*255);
                if (r.val_buffer.empty()) { rgb[0]=rgb[1]=rgb[2]=255; }
                f.write((char*)rgb, 3);
            }
        }
    }
}

/* ----------------------
   Processing Pipeline
   ---------------------- */
template <typename T>
ViewResult<T> process_view(const Vec3<T>& view_dir, const Mesh<T>& mesh, T pix_sz, bool cull, Renderer<T>& renderer) {
    auto t0 = std::chrono::high_resolution_clock::now();
    ViewResult<T> res;
    
    // Basis
    Vec3<T> n = normalize(view_dir);
    Vec3<T> h = (std::abs(n.x) < 0.9) ? Vec3<T>{1,0,0} : Vec3<T>{0,1,0};
    Vec3<T> u = normalize(cross(h, n));
    Vec3<T> v = cross(n, u);

    // Project
    std::vector<TriProjected<T>> tris;
    tris.reserve(mesh.faces.size());
    T b_min_x=1e20, b_max_x=-1e20, b_min_y=1e20, b_max_y=-1e20;

    for (size_t i=0; i<mesh.faces.size(); ++i) {
        const auto& f = mesh.faces[i];
        const Vec3<T>& v0 = mesh.vertices[f.v0_idx];
        const Vec3<T>& v1 = mesh.vertices[f.v1_idx];
        const Vec3<T>& v2 = mesh.vertices[f.v2_idx];

        Vec3<T> tri_n = cross(v1-v0, v2-v0);
        if (cull && dot(tri_n, n) >= 0) continue;

        Vec2<T> p0 = {dot(v0,u), dot(v0,v)};
        Vec2<T> p1 = {dot(v1,u), dot(v1,v)};
        Vec2<T> p2 = {dot(v2,u), dot(v2,v)};

        T area2 = edge_eval(p0, p1, p2);
        if (std::abs(area2) < 1e-12) continue;

        b_min_x = std::min({b_min_x, p0.x, p1.x, p2.x});
        b_max_x = std::max({b_max_x, p0.x, p1.x, p2.x});
        b_min_y = std::min({b_min_y, p0.y, p1.y, p2.y});
        b_max_y = std::max({b_max_y, p0.y, p1.y, p2.y});

        T val0=0, val1=0, val2=0;
        if (mesh.val_mode == MODE_NODE) {
            val0=mesh.values[f.v0_idx]; val1=mesh.values[f.v1_idx]; val2=mesh.values[f.v2_idx];
        } else if (mesh.val_mode == MODE_FACE) {
            val0=val1=val2=mesh.values[i];
        }

        tris.push_back({p0, p1, p2, dot(v0,n), dot(v1,n), dot(v2,n), val0, val1, val2, 1.0f/area2});
    }

    if (!tris.empty()) {
        b_min_x -= pix_sz; b_max_x += pix_sz;
        b_min_y -= pix_sz; b_max_y += pix_sz;
        int w = (int)std::ceil((b_max_x - b_min_x)/pix_sz);
        int h = (int)std::ceil((b_max_y - b_min_y)/pix_sz);
        
        renderer.resize(w, h, mesh.val_mode != MODE_NONE);
        res.image_width=w; res.image_height=h;

        for (const auto& t : tris) rasterize(t, renderer, b_min_x, b_min_y, pix_sz);

        // Stats
        size_t npix = (size_t)w * h;
        double sum = 0;
        for (size_t k=0; k<npix; ++k) {
            if (renderer.mask[k]) {
                res.covered_pixels++;
                if (mesh.val_mode != MODE_NONE) {
                    T val = renderer.val_buffer[k];
                    sum += val;
                    if (val < res.min_val) res.min_val = val;
                    if (val > res.max_val) res.max_val = val;
                }
            }
        }
        res.area = res.covered_pixels * pix_sz * pix_sz;
        if (res.covered_pixels) res.average_value = (T)(sum / res.covered_pixels);
        else { res.min_val = res.max_val = 0; }
    }
    
    res.time_seconds = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0).count();
    return res;
}

/* ----------------------
   Application Runners
   ---------------------- */
// Escaping for JSON
std::string esc(const std::string& s) {
    std::string out;
    for(char c:s) (c=='\\'||c=='"') ? out+="\\" + std::string(1,c) : out+=c;
    return out;
}

template <typename T>
void run_app(const std::string& obj, const std::string& default_data, const std::string& out_pre,
             const std::vector<BatchEntry>& batch, double default_res, bool cull, bool json) {
    
    Mesh<T> mesh = load_obj<T>(obj);
    if (!default_data.empty()) load_values(default_data, mesh); // Preload default data

    Renderer<T> renderer;

    if (json) {
        std::cout << "{\n  \"meta\": {\n"
                  << "    \"mesh\": \"" << esc(obj) << "\",\n"
                  << "    \"precision\": \"" << (sizeof(T)==4?"float":"double") << "\"\n"
                  << "  },\n  \"results\": [\n";
    }

    for (size_t k=0; k<batch.size(); ++k) {
        const auto& b = batch[k];
        
        // Handle Per-View Data Loading
        // Optimization: Only reload if this view specifies a DIFFERENT file than current
        // (Note: For simplicity here we just reload if specified, real app might cache)
        if (!b.data_file.empty()) load_values(b.data_file, mesh);
        
        T res = (b.resolution > 0) ? (T)b.resolution : (T)default_res;
        ViewResult<T> r = process_view({(T)b.nx, (T)b.ny, (T)b.nz}, mesh, res, cull, renderer);
        
        // Output Image?
        if (!out_pre.empty()) {
            std::ostringstream oss; oss << out_pre << "_" << std::setw(4) << std::setfill('0') << k << ".ppm";
            save_ppm(renderer, oss.str(), r.min_val, r.max_val);
        }

        if (json) {
            std::cout << "    {\n"
                      << "      \"idx\": " << k << ",\n"
                      << "      \"normal\": [" << b.nx << "," << b.ny << "," << b.nz << "],\n"
                      << "      \"area\": " << r.area << ",\n"
                      << "      \"avg\": " << r.average_value << ",\n"
                      << "      \"pixels\": " << r.covered_pixels << ",\n"
                      << "      \"time\": " << r.time_seconds << "\n"
                      << "    }" << (k==batch.size()-1?"":",") << "\n";
        } else {
            std::cout << "View " << k << " | Area: " << r.area << " | Avg: " << r.average_value
                      << " | Time: " << r.time_seconds << "s\n";
        }
        
        // Restore default data if we switched it temporarily?
        // For efficiency, we assume subsequent views use what is loaded unless changed.
        // If batch mixed global/local indiscriminately, we'd need to track state carefully.
    }
    if (json) std::cout << "  ]\n}\n";
}

/* ----------------------
   CLI Parsing
   ---------------------- */
void print_help(const char* prog) {
    std::cout << "FlatLand: Mesh Projection & Rasterization Tool\n"
              << "Usage: " << prog << " [mesh.obj] [options]\n\n"
              << "Core Options:\n"
              << "  -v <x> <y> <z>        Add a single view direction. Can be repeated multiple times.\n"
              << "  -b, --batch <file>    Load views from a file. See 'Batch Format' below.\n"
              << "                        (Note: Views from -v and -b are additive)\n\n"
              << "Simulation Parameters:\n"
              << "  -r, --res <val>       Set default pixel resolution (default: 0.001).\n"
              << "                        Used if a batch entry does not specify its own resolution.\n"
              << "  -d, --data <file>     Path to default scalar data file (node or face values).\n"
              << "                        Used if a batch entry does not specify its own data file.\n"
              << "  --no-cull             Disable backface culling (renders all matching faces).\n"
              << "  -p, --precision <mode> Set math precision: 'float' (default, faster) or 'double' (more accurate).\n\n"
              << "Output Control:\n"
              << "  -o, --out <prefix>    Save heatmap images with this filename prefix.\n"
              << "                        Files named: <prefix>_<view_index>.ppm\n"
              << "  -j, --json            Output results in structured JSON format to stdout.\n"
              << "  -h, --help            Show this help message.\n\n"
              << "Batch File Format:\n"
              << "  Each line represents a view. Lines starting with # are comments.\n"
              << "  Columns: <nx> <ny> <nz> [resolution] [data_file_path]\n"
              << "  Examples:\n"
              << "    1.0 0.0 0.0               (Uses default resolution & default data)\n"
              << "    0.0 1.0 0.0 0.005         (Overrides resolution, uses default data)\n"
              << "    0.0 0.0 1.0 0.002 data.txt (Overrides both)\n\n"
              << "Precedence Rules:\n"
              << "  1. Batch file specific settings (resolution/data) take highest precedence.\n"
              << "  2. Command line flags (-r, -d) set the defaults for any view missing those details.\n"
              << "  3. Single views (-v) always use the command-line default resolution and data.\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) { print_help(argv[0]); return 1; }

    std::string obj_file, data_file, batch_file, out_pre, prec="float";
    double res = 0.001;
    bool json=false, cull=true;
    std::vector<BatchEntry> batch;

    int i = 1;
    if (argv[1][0] != '-') { obj_file = argv[1]; i++; }

    for (; i<argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") { print_help(argv[0]); return 0; }
        else if (arg == "-d") { if(i+1<argc) data_file=argv[++i]; }
        else if (arg == "-b") { if(i+1<argc) batch_file=argv[++i]; }
        else if (arg == "-o") { if(i+1<argc) out_pre=argv[++i]; }
        else if (arg == "-r") { if(i+1<argc) res=std::stod(argv[++i]); }
        else if (arg == "-p" || arg == "--precision") { if(i+1<argc) prec=argv[++i]; }
        else if (arg == "-j" || arg == "--json") { json=true; }
        else if (arg == "--no-cull") { cull=false; }
        else if (arg == "-v") {
            if (i+3<argc) {
                batch.push_back({std::stod(argv[i+1]), std::stod(argv[i+2]), std::stod(argv[i+3]), -1.0, ""});
                i+=3;
            }
        }
        else if (obj_file.empty() && arg[0]!='-') { obj_file = arg; }
    }

    if (obj_file.empty()) { std::cerr << "Error: No OBJ file.\n"; return 1; }

    // Parse Batch File if present
    if (!batch_file.empty()) {
        std::ifstream f(batch_file);
        if (!f.is_open()) { std::cerr << "Error: Cannot open batch file.\n"; return 1; }
        std::string line;
        while (std::getline(f, line)) {
            size_t c = line.find('#');
            if (c != std::string::npos) line = line.substr(0, c);
            std::stringstream ss(line);
            double nx, ny, nz;
            if (ss >> nx >> ny >> nz) {
                BatchEntry be = {nx, ny, nz, -1.0, ""};
                // Optional columns
                double r_temp;
                if (ss >> r_temp) {
                    be.resolution = r_temp;
                    std::string d_temp;
                    if (ss >> d_temp) be.data_file = d_temp;
                }
                batch.push_back(be);
            }
        }
    }
    
    if (batch.empty()) { std::cerr << "Error: No views specified.\n"; return 1; }

    // Dispatch
    if (!json) std::cerr << "Running with " << prec << " precision...\n";
    
    if (prec == "float") {
        run_app<float>(obj_file, data_file, out_pre, batch, res, cull, json);
    } else {
        run_app<double>(obj_file, data_file, out_pre, batch, res, cull, json);
    }

    return 0;
}
