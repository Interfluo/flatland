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

#include "fl_batch.hpp"
#include "fl_mesh.hpp"
#include "fl_vec.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#include <clocale>
#include <cstddef>
#include <cstdlib>

namespace flatland {

/* ----------------------
   CLI Parsing
   ---------------------- */
void print_help(const char* prog) {
    std::cout << "FlatLand: Mesh Projection & Rasterization Tool\n"
              << "Usage: " << prog << " <mesh.obj|.stl> [options]\n\n"
              << "Meshes:\n"
              << "  OBJ and STL (binary or ASCII) are accepted, auto-detected by extension.\n\n"
              << "Core Options:\n"
              << "  -v, --view <x> <y> <z>   Add a view direction (vector). Repeatable.\n"
              << "  -a, --angle <az> <el>    Add a view by azimuth/elevation in degrees.\n"
              << "                           az sweeps around +Z from +X; el rises from XY.\n"
              << "  -b, --batch <file>       Load views from a file (see 'Batch Format').\n"
              << "                           Views from -v/-a and -b are additive.\n\n"
              << "Analysis Parameters:\n"
              << "  -r, --res <val>          Default pixel resolution (default: 0.001).\n"
              << "                           Used when a batch entry omits its own resolution.\n"
              << "  -d, --data <file[@col]>  Default scalar field file (see 'Field Files').\n"
              << "      --field-mode <m>     Field is 'node', 'face', or 'auto' (default auto,\n"
              << "                           inferred from row count).\n"
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
              << "Field Files:\n"
              << "  One row per mesh entity (vertex or face); one or more columns. Each\n"
              << "  column is a timestep, so a whole time series fits in ONE file. Select a\n"
              << "  column with a trailing '@<col>' (default 0), e.g. 'fields.txt@7'.\n"
              << "  A plain one-value-per-line file is just the single-column case.\n\n"
              << "Batch File Format:\n"
              << "  One view per line; '#' starts a comment. Two line forms:\n"
              << "    <nx> <ny> <nz> [resolution] [data[@col]]   (direction vector)\n"
              << "    a <az> <el>    [resolution] [data[@col]]   (azimuth/elevation)\n"
              << "  resolution and data are optional and order-independent, e.g.:\n"
              << "    1 0 0                      (default resolution & data)\n"
              << "    a 45 30 0.005              (angle view, custom resolution)\n"
              << "    0 1 0 fields.txt@12        (default resolution, column 12 of a series)\n";
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

} // namespace flatland

using namespace flatland;

int main(int argc, char* argv[]) {
    // Force the C locale so numeric parsing (strtod/sscanf) always uses '.' as the
    // decimal separator, regardless of the user's environment locale.
    std::setlocale(LC_ALL, "C");

    if (argc < 2) { print_help(argv[0]); return 1; }

    std::string obj_file, data_file, batch_file, out_pre, prec="float";
    double res = 0.001;
    bool json=false, cull=true;
    unsigned threads = 0; // 0 => auto
    ValueMode forced_mode = MODE_NONE; // --field-mode auto by default
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
            else if (a == "--field-mode") {
                std::string m = need_value(i,argc,argv,a);
                if      (m == "node") forced_mode = MODE_NODE;
                else if (m == "face") forced_mode = MODE_FACE;
                else if (m == "auto") forced_mode = MODE_NONE;
                else throw std::runtime_error("--field-mode must be 'node', 'face', or 'auto', got '" + m + "'");
            }
            else if (a == "-v" || a == "--view") {
                if (i+3 >= argc) throw std::runtime_error("option '" + a + "' requires three numbers: x y z");
                double x = parse_double(argv[i+1], "-v x");
                double y = parse_double(argv[i+2], "-v y");
                double z = parse_double(argv[i+3], "-v z");
                if (x==0 && y==0 && z==0) throw std::runtime_error("view direction cannot be the zero vector");
                batch.push_back({x, y, z, -1.0, ""});
                i += 3;
            }
            else if (a == "-a" || a == "--angle") {
                if (i+2 >= argc) throw std::runtime_error("option '" + a + "' requires two numbers: azimuth elevation (degrees)");
                double az = parse_double(argv[i+1], "-a azimuth");
                double el = parse_double(argv[i+2], "-a elevation");
                Vec3<double> d = angle_to_dir(az, el);
                batch.push_back({d.x, d.y, d.z, -1.0, ""});
                i += 2;
            }
            else if (!a.empty() && a[0] == '-') {
                throw std::runtime_error("unknown option '" + a + "' (try --help)");
            }
            else if (obj_file.empty()) { obj_file = a; }
            else throw std::runtime_error("unexpected argument '" + a + "' (mesh already set to '" + obj_file + "')");
        }

        if (obj_file.empty()) throw std::runtime_error("no mesh file specified (.obj or .stl)");
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
                std::string first;
                if (!(ss >> first)) continue;   // blank / comment-only line

                // A line is either a direction vector "<nx> <ny> <nz> ..." or an
                // angle "a <az> <el> ..." (degrees). Detect by the first token.
                double nx, ny, nz;
                if (first == "a" || first == "A" || first == "angle") {
                    double az, el;
                    if (!(ss >> az >> el))
                        throw std::runtime_error("batch line " + std::to_string(lineno) + ": angle view expects 'a <azimuth> <elevation> ...'");
                    Vec3<double> d = angle_to_dir(az, el);
                    nx = d.x; ny = d.y; nz = d.z;
                } else {
                    char* endp; nx = std::strtod(first.c_str(), &endp);
                    if (*endp != '\0' || !(ss >> ny >> nz))
                        throw std::runtime_error("batch line " + std::to_string(lineno) + ": expected '<nx> <ny> <nz> ...' or 'a <az> <el> ...'");
                    if (nx==0 && ny==0 && nz==0)
                        throw std::runtime_error("batch line " + std::to_string(lineno) + ": zero view vector");
                }

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
            }
        }

        if (batch.empty()) throw std::runtime_error("no views specified (use -v or -b)");

        unsigned hw = std::thread::hardware_concurrency();
        if (threads == 0) threads = hw ? hw : 1;

        if (!json) std::cerr << "Running " << batch.size() << " view(s) with " << prec
                             << " precision on " << threads << " thread(s)...\n";

        if (prec == "float") run_app<float>(obj_file, data_file, out_pre, batch, res, cull, json, threads, forced_mode);
        else                 run_app<double>(obj_file, data_file, out_pre, batch, res, cull, json, threads, forced_mode);

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
