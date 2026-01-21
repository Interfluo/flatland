//
//  main.cpp
//  FlatLand
//
//  Created by James Crouse on 12/14/25.
//

#include <iostream>
#include <vector>
#include <cstdint>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <limits>

/* ----------------------
   Custom Types
   ---------------------- */
struct Vec3 { double x = 0.0, y = 0.0, z = 0.0; };
struct Vec2 { double x = 0.0, y = 0.0; };

struct Triangle3D { Vec3 v0, v1, v2; };
struct Triangle2D {
    Vec2 v0, v1, v2;
    double d0, d1, d2;  // depths
};

struct Mesh { std::vector<Triangle3D> triangles; };

struct ViewFrame {
    Vec3 u;  // x-axis
    Vec3 v;  // y-axis
    Vec3 n;  // normalized view normal
};

struct Bounds2D {
    double xmin =  std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    double ymin =  std::numeric_limits<double>::max();
    double ymax = -std::numeric_limits<double>::max();
};

struct Image {
    int Nx = 0, Ny = 0;
    double xmin = 0.0, ymin = 0.0, pixel_size = 0.0;
    std::vector<double> zbuffer;
    std::vector<uint8_t> mask;
};

struct ViewResult {
    double area = 0.0;
    double time_seconds = 0.0;
    int covered_pixels = 0;
    int visible_triangles = 0;
    int image_width = 0;
    int image_height = 0;
    std::string output_image;  // empty if no image written
};

/* ----------------------
   Vector Math
   ---------------------- */
inline Vec3 operator+(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 operator-(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec2 operator-(Vec2 a, Vec2 b) { return {a.x - b.x, a.y - b.y}; }
inline Vec3 operator*(double s, Vec3 v) { return {s * v.x, s * v.y, s * v.z}; }
inline double dot(Vec3 a, Vec3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline double dot(Vec2 a, Vec2 b) { return a.x*b.x + a.y*b.y; }
inline Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}
inline double norm(Vec3 v) { return std::sqrt(dot(v, v)); }
inline Vec3 normalize(Vec3 v) { double n = norm(v); return n > 0.0 ? (1.0/n)*v : Vec3{}; }
Vec3 triangle_normal(const Triangle3D& t) { return cross(t.v1 - t.v0, t.v2 - t.v0); }

/* ----------------------
   Mesh IO
   ---------------------- */
Mesh load_obj(const std::string& filename) {
    Mesh mesh;
    std::vector<Vec3> vertices;
    std::ifstream file(filename);
    if (!file.is_open()) throw std::invalid_argument("Error: Could not open OBJ file '" + filename + "'");

    std::string line;
    size_t line_number = 0;
    while (std::getline(file, line)) {
        ++line_number;
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;

        if (prefix == "v") {
            Vec3 v;
            if (!(iss >> v.x >> v.y >> v.z))
                std::cerr << "Warning: Malformed vertex on line " << line_number << "\n";
            else
                vertices.push_back(v);
        } else if (prefix == "f") {
            std::vector<int> face;
            std::string token;
            while (iss >> token) {
                size_t slash = token.find('/');
                std::string idx_str = (slash == std::string::npos) ? token : token.substr(0, slash);
                if (idx_str.empty()) continue;
                int idx = std::stoi(idx_str);
                if (idx > 0) idx -= 1;
                else idx += static_cast<int>(vertices.size());
                if (idx < 0 || static_cast<size_t>(idx) >= vertices.size()) {
                    std::cerr << "Warning: Invalid vertex index on line " << line_number << "\n";
                    face.clear();
                    break;
                }
                face.push_back(idx);
            }
            if (face.size() == 3) {
                mesh.triangles.push_back({vertices[face[0]], vertices[face[1]], vertices[face[2]]});
            } else if (face.size() == 4) {
                mesh.triangles.push_back({vertices[face[0]], vertices[face[1]], vertices[face[2]]});
                mesh.triangles.push_back({vertices[face[0]], vertices[face[2]], vertices[face[3]]});
            } else if (!face.empty()) {
                std::cerr << "Warning: Skipping face with " << face.size()
                          << " vertices on line " << line_number << " (only tris/quads supported)\n";
            }
        }
    }
    if (vertices.empty()) throw std::invalid_argument("Error: No vertices loaded.");
    return mesh;
}

/* ----------------------
   Rasterization Utilities
   ---------------------- */
bool point_in_triangle(const Vec2& pt, const Vec2& v0, const Vec2& v1, const Vec2& v2,
                       double& alpha, double& beta, double& gamma) {
    Vec2 e0 = v1 - v0, e1 = v2 - v0, ep = pt - v0;
    double d00 = dot(e0, e0), d01 = dot(e0, e1), d11 = dot(e1, e1);
    double d0p = dot(e0, ep), d1p = dot(e1, ep);
    double denom = d00*d11 - d01*d01;
    if (std::abs(denom) < 1e-20) return false;
    double inv = 1.0 / denom;
    beta  = (d11*d0p - d01*d1p) * inv;
    gamma = (d00*d1p - d01*d0p) * inv;
    alpha = 1.0 - beta - gamma;
    const double eps = 1e-10;
    return alpha > -eps && beta > -eps && gamma > -eps;
}

size_t pixel_index(int ix, int iy, const Image& img) { return static_cast<size_t>(iy) * img.Nx + ix; }
Vec2 pixel_center(int ix, int iy, const Image& img) {
    return {img.xmin + (ix + 0.5)*img.pixel_size,
            img.ymin + (iy + 0.5)*img.pixel_size};
}

void rasterize_triangle(const Triangle2D& tri, Image& img) {
    double min_x = std::min({tri.v0.x, tri.v1.x, tri.v2.x});
    double max_x = std::max({tri.v0.x, tri.v1.x, tri.v2.x});
    double min_y = std::min({tri.v0.y, tri.v1.y, tri.v2.y});
    double max_y = std::max({tri.v0.y, tri.v1.y, tri.v2.y});

    int ix_min = std::max(0, static_cast<int>(std::floor((min_x - img.xmin)/img.pixel_size)) - 1);
    int ix_max = std::min(img.Nx - 1, static_cast<int>(std::floor((max_x - img.xmin)/img.pixel_size)) + 1);
    int iy_min = std::max(0, static_cast<int>(std::floor((min_y - img.ymin)/img.pixel_size)) - 1);
    int iy_max = std::min(img.Ny - 1, static_cast<int>(std::floor((max_y - img.ymin)/img.pixel_size)) + 1);

    for (int iy = iy_min; iy <= iy_max; ++iy) {
        for (int ix = ix_min; ix <= ix_max; ++ix) {
            Vec2 pt = pixel_center(ix, iy, img);
            double alpha, beta, gamma;
            if (point_in_triangle(pt, tri.v0, tri.v1, tri.v2, alpha, beta, gamma)) {
                double depth = alpha*tri.d0 + beta*tri.d1 + gamma*tri.d2;
                size_t idx = pixel_index(ix, iy, img);
                if (depth < img.zbuffer[idx]) {
                    img.zbuffer[idx] = depth;
                    img.mask[idx] = 1;
                }
            }
        }
    }
}

/* ----------------------
   Misc Utilities
   ---------------------- */
ViewFrame make_view_frame(const Vec3& view_normal) {
    ViewFrame f;
    f.n = normalize(view_normal);
    Vec3 helper = (std::abs(f.n.x) < 0.9) ? Vec3{1.0,0.0,0.0} : Vec3{0.0,1.0,0.0};
    f.u = normalize(cross(helper, f.n));
    f.v = cross(f.n, f.u);
    return f;
}

void write_ppm_silhouette(const Image& img, const std::string& filename) {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open '" << filename << "' for writing PPM\n";
        return;
    }
    file << "P6\n" << img.Nx << " " << img.Ny << "\n255\n";
    for (int iy = img.Ny-1; iy >= 0; --iy) {
        for (int ix = 0; ix < img.Nx; ++ix) {
            uint8_t v = img.mask[pixel_index(ix, iy, img)] ? 0 : 255;
            uint8_t pixel[3] = {v, v, v};
            file.write(reinterpret_cast<const char*>(pixel), 3);
        }
    }
}

std::string get_output_default(const std::string& input_path) {
    size_t last_slash = input_path.find_last_of("\\/");
    std::string name = (last_slash == std::string::npos) ? input_path : input_path.substr(last_slash + 1);
    size_t dot = name.rfind('.');
    if (dot != std::string::npos) name = name.substr(0, dot);
    return name + ".ppm";
}

Bounds2D make_empty_bounds() {
    Bounds2D b;
    b.xmin = b.ymin =  std::numeric_limits<double>::max();
    b.xmax = b.ymax = -std::numeric_limits<double>::max();
    return b;
}

void expand_bounds(Bounds2D& b, const Vec2& p) {
    b.xmin = std::min(b.xmin, p.x);
    b.xmax = std::max(b.xmax, p.x);
    b.ymin = std::min(b.ymin, p.y);
    b.ymax = std::max(b.ymax, p.y);
}

/* ----------------------
   Core Per-View Processing (modularized)
   ---------------------- */
ViewResult process_single_view(const Vec3& view_normal,
                               const Mesh& mesh,
                               double pixel_size,
                               bool back_face_cull,
                               bool need_image,
                               const std::string& image_prefix,
                               size_t view_index,
                               bool is_batch) {
    ViewResult result;

    auto start = std::chrono::high_resolution_clock::now();

    ViewFrame frame = make_view_frame(view_normal);

    Bounds2D bounds = make_empty_bounds();
    std::vector<Triangle2D> projected;

    auto project = [&](const Vec3& p) { return Vec2{dot(p, frame.u), dot(p, frame.v)}; };
    auto depth   = [&](const Vec3& p) { return dot(p, frame.n); };

    for (const auto& tri : mesh.triangles) {
        if (back_face_cull && dot(triangle_normal(tri), frame.n) >= 0.0) continue;

        Vec2 p0 = project(tri.v0), p1 = project(tri.v1), p2 = project(tri.v2);
        double d0 = depth(tri.v0), d1 = depth(tri.v1), d2 = depth(tri.v2);

        projected.push_back({p0, p1, p2, d0, d1, d2});
        expand_bounds(bounds, p0);
        expand_bounds(bounds, p1);
        expand_bounds(bounds, p2);
        ++result.visible_triangles;
    }

    bool has_geometry = (bounds.xmax > bounds.xmin && bounds.ymax > bounds.ymin);
    if (has_geometry) {
        bounds.xmin -= pixel_size; bounds.xmax += pixel_size;
        bounds.ymin -= pixel_size; bounds.ymax += pixel_size;

        double width  = bounds.xmax - bounds.xmin;
        double height = bounds.ymax - bounds.ymin;

        result.image_width  = static_cast<int>(std::ceil(width / pixel_size));
        result.image_height = static_cast<int>(std::ceil(height / pixel_size));

        if (result.image_width > 0 && result.image_height > 0) {
            Image img;
            img.Nx = result.image_width;
            img.Ny = result.image_height;
            img.xmin = bounds.xmin;
            img.ymin = bounds.ymin;
            img.pixel_size = pixel_size;
            img.zbuffer.assign(static_cast<size_t>(img.Nx) * img.Ny, std::numeric_limits<double>::max());
            img.mask.assign(static_cast<size_t>(img.Nx) * img.Ny, 0);

            for (const auto& t : projected) rasterize_triangle(t, img);

            result.covered_pixels = 0;
            for (uint8_t m : img.mask) result.covered_pixels += m;
            result.area = result.covered_pixels * pixel_size * pixel_size;

            if (need_image) {
                std::ostringstream oss;
                oss << image_prefix;
                if (is_batch) oss << "_" << std::setw(4) << std::setfill('0') << view_index;
                oss << ".ppm";
                result.output_image = oss.str();
                write_ppm_silhouette(img, result.output_image);
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    result.time_seconds = std::chrono::duration<double>(end - start).count();

    return result;
}

/* ----------------------
   Help & Main
   ---------------------- */
void print_help(const char* prog) {
    std::cout << "Usage: " << prog << " <input.obj> [options]\n\n"
              << "Computes visible projected area of an OBJ mesh from one or more view directions.\n\n"
              << "Options:\n"
              << "  -h, --help                Show this help\n"
              << "  -v, --verbose             Verbose output\n"
              << "  --view <x> <y> <z>        Single view normal (default: 1 0 0)\n"
              << "  --views-file <file>       Text file with multiple view normals (x y z per line)\n"
              << "  -p, --pixel-size <val>    Pixel size (default: 0.001)\n"
              << "  -o, --output <prefix>     Output PPM prefix (batch: prefix_0000.ppm etc.)\n"
              << "  --no-image                Disable PPM output\n"
              << "  --json                    JSON output\n"
              << "  --no-bfc                  Disable back-face culling\n";
}

int main(int argc, const char* argv[]) {
    if (argc < 2) { std::cerr << "Error: Missing input OBJ file.\n"; print_help(argv[0]); return 1; }

    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "-h" || std::string(argv[i]) == "--help") { print_help(argv[0]); return 0; }

    std::string mesh_filename = argv[1];
    std::string views_file, output_prefix;
    Vec3 single_view{1.0, 0.0, 0.0};
    bool single_set = false;
    double pixel_size = 0.001;
    bool verbose = false, write_image = true, json = false, bfc = true;

    for (int i = 2; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "-v" || arg == "--verbose") verbose = true;
        else if (arg == "--view") {
            if (i+3 >= argc) { std::cerr << "Error: --view needs 3 values\n"; return 1; }
            single_view = {std::stod(argv[++i]), std::stod(argv[++i]), std::stod(argv[++i])};
            single_set = true;
        } else if (arg == "--views-file") {
            if (i+1 >= argc) { std::cerr << "Error: --views-file needs filename\n"; return 1; }
            views_file = argv[++i];
        } else if (arg == "-p" || arg == "--pixel-size") {
            if (i+1 >= argc) { std::cerr << "Error: --pixel-size needs value\n"; return 1; }
            pixel_size = std::stod(argv[++i]);
        } else if (arg == "-o" || arg == "--output") {
            if (i+1 >= argc) { std::cerr << "Error: --output needs filename\n"; return 1; }
            output_prefix = argv[++i];
            write_image = true;
        } else if (arg == "--no-image") write_image = false;
        else if (arg == "--json") json = true;
        else if (arg == "--no-bfc") bfc = false;
        else { std::cerr << "Unknown option: " << arg << "\n"; return 1; }
    }

    // Load views
    std::vector<Vec3> views;
    if (!views_file.empty()) {
        std::ifstream f(views_file);
        if (!f) { std::cerr << "Error: Cannot open views file\n"; return 1; }
        std::string line; size_t ln = 0;
        while (std::getline(f, line)) {
            ++ln;
            if (line.empty() || line[0] == '#') continue;
            std::istringstream iss(line);
            Vec3 v;
            if (iss >> v.x >> v.y >> v.z && norm(v) > 0.0) views.push_back(v);
            else std::cerr << "Warning: Invalid view on line " << ln << "\n";
        }
        if (views.empty()) { std::cerr << "Error: No valid views in file\n"; return 1; }
    } else {
        Vec3 v = single_set ? single_view : Vec3{1.0,0.0,0.0};
        if (norm(v) == 0.0) { std::cerr << "Error: Zero view normal\n"; return 1; }
        views.push_back(v);
    }

    bool batch = views.size() > 1;
    bool need_image = write_image && !(batch && output_prefix.empty());
    if (need_image && output_prefix.empty())
        output_prefix = get_output_default(mesh_filename);
    size_t dot = output_prefix.rfind('.');
    if (dot != std::string::npos && output_prefix.substr(dot) == ".ppm")
        output_prefix = output_prefix.substr(0, dot);

    Mesh mesh = load_obj(mesh_filename);
    if (verbose && !json) std::cout << "Loaded mesh with " << mesh.triangles.size() << " triangles\n";

    if (json && batch) std::cout << "[\n";
    bool first = true;
    double total_time = 0.0;

    for (size_t i = 0; i < views.size(); ++i) {
        ViewResult res = process_single_view(views[i], mesh, pixel_size, bfc,
                                             need_image, output_prefix, i, batch);

        total_time += res.time_seconds;

        if (!json) {
            if (verbose) {
                std::cout << "View " << (i+1) << "/" << views.size()
                          << " (" << views[i].x << " " << views[i].y << " " << views[i].z << ")\n";
                std::cout << "  Visible triangles: " << res.visible_triangles << " / " << mesh.triangles.size() << "\n";
                std::cout << "  Resolution: " << res.image_width << "x" << res.image_height << "\n";
                std::cout << "  Covered pixels: " << res.covered_pixels << "\n";
                std::cout << "  Time: " << std::fixed << std::setprecision(3) << res.time_seconds << " s\n";
            }
            std::cout << "--> Projected area: " << std::fixed << std::setprecision(12)
                      << res.area << " units²\n";
            if (!res.output_image.empty())
                std::cout << "Wrote silhouette to '" << res.output_image << "'\n";
            if (batch) std::cout << "\n";
        } else {
            if (batch && !first) std::cout << ",\n";
            first = false;
            std::cout << "  {\n"
                      << "    \"view_normal\": [" << views[i].x << "," << views[i].y << "," << views[i].z << "],\n"
                      << "    \"area\": " << std::fixed << std::setprecision(12) << res.area << ",\n"
                      << "    \"time_seconds\": " << std::setprecision(6) << res.time_seconds << ",\n"
                      << "    \"covered_pixels\": " << res.covered_pixels << ",\n"
                      << "    \"visible_triangles\": " << res.visible_triangles << ",\n"
                      << "    \"total_triangles\": " << mesh.triangles.size() << ",\n"
                      << "    \"image_width\": " << res.image_width << ",\n"
                      << "    \"image_height\": " << res.image_height << ",\n"
                      << "    \"pixel_size\": " << pixel_size;
            if (!res.output_image.empty())
                std::cout << ",\n    \"output_image\": \"" << res.output_image << "\"";
            std::cout << "\n  }";
        }
    }

    if (json && batch) std::cout << "\n]\n";
    if (!json && batch)
        std::cout << "Total time for " << views.size() << " views: "
                  << std::fixed << std::setprecision(3) << total_time << " s\n";

    return 0;
}
