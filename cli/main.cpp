// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Interfluo
//
// FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
// for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

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
//   - PNG heatmap export, raw .npy value export, and structured JSON output
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
#include <cctype>
#include <cmath>
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
              << "  -o, --out <prefix>       Save PNG heatmaps as <prefix>_<idx>.png.\n"
              << "      --npy <prefix>       Save raw per-pixel field values as\n"
              << "                           <prefix>_<idx>.npy (NumPy float64, NaN off the\n"
              << "                           silhouette). Needs a field; row 0 is the BOTTOM\n"
              << "                           row, so plot with origin='lower'.\n"
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

// Upper bound on -t. The clamp against the view count still applies; this exists
// so that an unbounded or negative request cannot reach std::thread at all.
// Previously `-t -1` wrapped to 4294967295 and `-t 100000` was taken literally,
// and on a large batch the thread constructor threw while a vector of joinable
// threads was still live, which terminates the process.
static const unsigned MAX_THREADS = 4096;

static unsigned parse_threads(const std::string& s) {
    const double d = parse_double(s, "-t/--threads");
    if (!std::isfinite(d) || d < 0 || d > (double)MAX_THREADS || d != std::floor(d))
        throw std::runtime_error("--threads must be a whole number from 0 to " +
                                 std::to_string(MAX_THREADS) + " (0 means one per core), got '" + s + "'");
    return (unsigned)d;
}

// POSIX '/', plus Windows '\' and "C:" so a batch written on either platform
// resolves the same way.
static bool is_absolute_path(const std::string& p) {
    if (p.empty()) return false;
    if (p[0] == '/' || p[0] == '\\') return true;
    return p.size() >= 2 && p[1] == ':' && std::isalpha((unsigned char)p[0]);
}

static std::string resolve_batch_path(const std::string& dir, const std::string& p) {
    return is_absolute_path(p) ? p : dir + p;
}

static bool readable(const std::string& p) {
    if (p.empty()) return false;
    std::ifstream f(p, std::ios::binary);
    return f.good();
}

// True when a batch token names a data file: either directly, or in the
// "<existing file>@<column>" form.
static bool names_data_file(const std::string& p) {
    if (readable(p)) return true;
    const size_t at = p.find_last_of('@');
    if (at == std::string::npos || at + 1 >= p.size()) return false;
    if (p.find_first_not_of("0123456789", at + 1) != std::string::npos) return false;
    return readable(p.substr(0, at));
}

} // namespace flatland

using namespace flatland;

int main(int argc, char* argv[]) {
    // Force the C locale so numeric parsing (strtod/sscanf) always uses '.' as the
    // decimal separator, regardless of the user's environment locale.
    std::setlocale(LC_ALL, "C");

    if (argc < 2) { print_help(argv[0]); return 1; }

    std::string obj_file, data_file, batch_file, out_pre, npy_pre, prec="float";
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
            else if (a == "--npy")                    { npy_pre    = need_value(i,argc,argv,a); }
            else if (a == "-r" || a == "--res")       { res = parse_double(need_value(i,argc,argv,a), "-r/--res"); }
            else if (a == "-p" || a == "--precision") { prec = need_value(i,argc,argv,a); }
            else if (a == "-t" || a == "--threads")   { threads = parse_threads(need_value(i,argc,argv,a)); }
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
                if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
                    throw std::runtime_error("view direction must be finite");
                batch.push_back({x, y, z, -1.0, ""});
                i += 3;
            }
            else if (a == "-a" || a == "--angle") {
                if (i+2 >= argc) throw std::runtime_error("option '" + a + "' requires two numbers: azimuth elevation (degrees)");
                double az = parse_double(argv[i+1], "-a azimuth");
                double el = parse_double(argv[i+2], "-a elevation");
                if (!std::isfinite(az) || !std::isfinite(el))
                    throw std::runtime_error("azimuth and elevation must be finite");
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
        // NaN and +inf both slip past a bare `res <= 0`, then surface far away
        // as an allocator error naming a std::vector internal.
        if (!std::isfinite(res) || res <= 0)
            throw std::runtime_error("resolution must be a positive, finite number, got '" +
                                     std::to_string(res) + "'");

        // Parse batch file (additive with any -v views)
        if (!batch_file.empty()) {
            // Binary for the same reason the mesh loader is: the tokenizer below
            // treats '\r' as whitespace, so reading the raw bytes makes a CRLF
            // batch file parse identically on Windows and on POSIX.
            std::ifstream f(batch_file, std::ios::binary);
            if (!f.is_open()) throw std::runtime_error("cannot open batch file '" + batch_file + "'");
            // Relative data-file paths in a batch are resolved against the batch file's
            // directory, so a committed time-series case works from any CWD.
            std::string batch_dir;
            size_t slash = batch_file.find_last_of("/\\");
            if (slash != std::string::npos) batch_dir = batch_file.substr(0, slash + 1);

            std::string line;
            int lineno = 0;
            while (std::getline(f, line)) {
                ++lineno;

                // Split into tokens, stopping at a token that STARTS with '#'.
                // Cutting the line at the first '#' anywhere would truncate a
                // legitimate path such as "run#3/field.txt".
                std::vector<std::string> tok;
                {
                    std::istringstream ss(line);
                    std::string t;
                    while (ss >> t) { if (t[0] == '#') break; tok.push_back(t); }
                }
                if (tok.empty()) continue;   // blank / comment-only line

                auto err = [&](const std::string& what) {
                    throw std::runtime_error("batch line " + std::to_string(lineno) + ": " + what);
                };

                // A line is either a direction vector "<nx> <ny> <nz> ..." or an
                // angle "a <az> <el> ..." (degrees). Detect by the first token.
                double nx, ny, nz;
                size_t ti;
                if (tok[0] == "a" || tok[0] == "A" || tok[0] == "angle") {
                    if (tok.size() < 3) err("angle view expects 'a <azimuth> <elevation> ...'");
                    double az = parse_double(tok[1], "batch azimuth");
                    double el = parse_double(tok[2], "batch elevation");
                    if (!std::isfinite(az) || !std::isfinite(el)) err("azimuth and elevation must be finite");
                    Vec3<double> d = angle_to_dir(az, el);
                    nx = d.x; ny = d.y; nz = d.z;
                    ti = 3;
                } else {
                    if (tok.size() < 3) err("expected '<nx> <ny> <nz> ...' or 'a <az> <el> ...'");
                    nx = parse_double(tok[0], "batch nx");
                    ny = parse_double(tok[1], "batch ny");
                    nz = parse_double(tok[2], "batch nz");
                    if (nx == 0 && ny == 0 && nz == 0) err("zero view vector");
                    if (!std::isfinite(nx) || !std::isfinite(ny) || !std::isfinite(nz))
                        err("view direction must be finite");
                    ti = 3;
                }

                BatchEntry be = {nx, ny, nz, -1.0, ""};

                // The optional trailing tokens are [resolution] and [data[@col]],
                // in EITHER order, as the help text and README have always said.
                //
                // The filesystem gets the first vote: a token naming a readable
                // file is data. Deciding on "looks like a number" first meant a
                // time-series file named for its timestep ("0100", "2024") was
                // silently read as a resolution and the field never loaded at
                // all — a wrong answer with a zero exit status. Ambiguity is
                // unavoidable here, so it resolves toward the loud failure.
                for (; ti < tok.size(); ++ti) {
                    const std::string resolved = resolve_batch_path(batch_dir, tok[ti]);
                    if (names_data_file(resolved)) {
                        if (!be.data_file.empty()) err("more than one data file given");
                        be.data_file = resolved;
                        continue;
                    }
                    char* endp = nullptr;
                    const char* cstr = tok[ti].c_str();
                    const double rv = std::strtod(cstr, &endp);
                    if (endp != cstr && *endp == '\0') {
                        if (!std::isfinite(rv) || rv <= 0)
                            err("resolution must be a positive, finite number, got '" + tok[ti] + "'");
                        if (be.resolution > 0) err("more than one resolution given");
                        be.resolution = rv;
                        continue;
                    }
                    // Neither an existing file nor a number. Treat it as a data
                    // path so the failure names the missing file, unless we
                    // already have one, in which case it is simply surplus.
                    if (!be.data_file.empty()) err("unexpected extra token '" + tok[ti] + "'");
                    be.data_file = resolved;
                }
                batch.push_back(be);
            }
        }
        if (batch.empty()) throw std::runtime_error("no views specified (use -v or -b)");

        unsigned hw = std::thread::hardware_concurrency();
        if (threads == 0) threads = hw ? hw : 1;

        if (prec == "float") run_app<float>(obj_file, data_file, out_pre, npy_pre, batch, res, cull, json, threads, forced_mode);
        else                 run_app<double>(obj_file, data_file, out_pre, npy_pre, batch, res, cull, json, threads, forced_mode);

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
