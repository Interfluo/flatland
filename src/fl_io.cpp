// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Interfluo
//
// FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
// for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

#include "fl_io.hpp"

#include "fl_locale.hpp"

#include <algorithm>
#include <fstream>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>

namespace flatland {

/* ----------------------
   IO & Loaders
   ----------------------
   Both loaders parse into double and hand a RawMesh to finalize<T>(), which
   validates, recenters and casts once. Parsing in double matters even for a
   float run: a vertex at 5e6 has already lost the mantissa bits that carry its
   detail by the time it is stored as a float, so recentering has to happen
   before the narrowing conversion, not after.

   The other rule enforced here is that malformed input is REPORTED. Silently
   skipping a vertex line shifts every subsequent 1-based face index by one and
   turns a parse failure into a plausible, wrong answer. */

namespace {

// A UTF-8 BOM is invisible in an editor but glued to the first token, so `v`
// arrives as "\xEF\xBB\xBFv" and the line is dropped.
void strip_bom(std::string& s) {
    if (s.size() >= 3 && (unsigned char)s[0] == 0xEF
                      && (unsigned char)s[1] == 0xBB
                      && (unsigned char)s[2] == 0xBF)
        s.erase(0, 3);
}

inline bool is_space(char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\n'; }

// First whitespace-delimited token, with the offset of whatever follows it.
// OBJ is free-form: leading whitespace is legal and plenty of exporters indent.
std::string first_token(const std::string& line, size_t& rest) {
    size_t b = 0;
    while (b < line.size() && is_space(line[b])) ++b;
    size_t e = b;
    while (e < line.size() && !is_space(line[e])) ++e;
    rest = e;
    return line.substr(b, e - b);
}

[[noreturn]] void fail(const std::string& file, size_t line, const std::string& what) {
    throw std::runtime_error("mesh '" + file + "' line " + std::to_string(line) + ": " + what);
}

RawMesh load_obj_raw(const std::string& filename) {
    // Binary, deliberately, even though OBJ is text. In text mode a Windows CRT
    // would eat the '\r' of a CRLF file and stop at an embedded 0x1A, so the
    // same file would parse into different bytes on different platforms. Every
    // reader below treats '\r' as whitespace (is_space), which is what the
    // CRLF case in tests/test_parsing.sh already exercises on POSIX, so reading
    // the raw bytes everywhere makes the two platforms agree by construction.
    // Numbers in an OBJ are written with '.', so strtod has to agree regardless
    // of the locale the host process is running in.
    CNumericScope c_numeric;

    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("cannot open mesh file '" + filename + "'");

    RawMesh mesh;
    std::string line;
    size_t lineno = 0;

    while (std::getline(file, line)) {
        ++lineno;
        if (lineno == 1) strip_bom(line);

        size_t rest = 0;
        const std::string kw = first_token(line, rest);
        if (kw.empty() || kw[0] == '#') continue;

        if (kw == "v") {
            // A geometric vertex is "v x y z [w]"; w is a rational weight we do
            // not use. Fewer than three coordinates is malformed, not ignorable.
            const char* p = line.c_str() + rest;
            double c[4];
            int got = 0;
            while (got < 4) {
                char* end;
                const double d = std::strtod(p, &end);
                if (end == p) break;
                c[got++] = d;
                p = end;
            }
            while (*p && is_space(*p)) ++p;
            if (got < 3)
                fail(filename, lineno, "vertex needs 3 coordinates, found " + std::to_string(got));
            if (*p)
                fail(filename, lineno, "unexpected text after the vertex coordinates");
            if (!std::isfinite(c[0]) || !std::isfinite(c[1]) || !std::isfinite(c[2]))
                fail(filename, lineno, "vertex coordinate is not a finite number");
            mesh.vertices.push_back({c[0], c[1], c[2]});

        } else if (kw == "f") {
            // Face entries are v, v/vt, v//vn or v/vt/vn; only the vertex index
            // matters. Negative indices are relative to the vertices seen so far.
            std::vector<int> idxs;
            const char* s = line.c_str() + rest;
            while (*s) {
                while (*s && is_space(*s)) ++s;
                if (!*s) break;
                char* end;
                // strtoll, not strtol: `long` is 32 bits on Windows, so a file
                // with a large index would saturate there and not on POSIX.
                const long long raw = std::strtoll(s, &end, 10);
                if (end == s)
                    fail(filename, lineno, "face vertex index is not a number");
                if (raw == 0)
                    fail(filename, lineno, "face vertex index 0 is invalid (OBJ indices are 1-based)");
                idxs.push_back(raw > 0 ? (int)(raw - 1)
                                       : (int)(raw + (long long)mesh.vertices.size()));
                while (*end && !is_space(*end)) ++end;   // skip /vt/vn
                s = end;
            }
            if (idxs.size() < 3)
                fail(filename, lineno, "face needs at least 3 vertices, found " + std::to_string(idxs.size()));
            for (size_t k = 2; k < idxs.size(); ++k)     // fan-triangulate
                mesh.faces.push_back({idxs[0], idxs[k-1], idxs[k]});
        }
        // Everything else (vn, vt, g, o, s, usemtl, mtllib, ...) is not geometry.
    }
    return mesh;
}

// STL, auto-detecting binary vs ASCII. STL stores independent per-triangle
// vertices (no shared indexing), so each triangle contributes three fresh ones.
RawMesh load_stl_raw(const std::string& filename) {
    CNumericScope c_numeric;          // ASCII STL parses with >>, also locale-sensitive
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("cannot open mesh file '" + filename + "'");
    const std::string buf((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

    RawMesh mesh;

    // A binary STL is exactly 84 + 50*ntri bytes. This is the robust
    // discriminator: the leading "solid" keyword is NOT reliable, since some
    // binary writers put it in the header.
    uint32_t ntri = 0;
    bool has_header = buf.size() >= 84;
    if (has_header) std::memcpy(&ntri, buf.data() + 80, 4);
    const bool binary = has_header && buf.size() == 84ull + 50ull * ntri;

    if (binary) {
        mesh.vertices.reserve((size_t)ntri * 3);
        mesh.faces.reserve(ntri);
        for (uint32_t t = 0; t < ntri; ++t) {
            const char* tri = buf.data() + 84 + (size_t)t * 50;
            float v[9];
            std::memcpy(v, tri + 12, 36);   // skip the 12-byte facet normal
            const int base = (int)mesh.vertices.size();
            mesh.vertices.push_back({v[0], v[1], v[2]});
            mesh.vertices.push_back({v[3], v[4], v[5]});
            mesh.vertices.push_back({v[6], v[7], v[8]});
            mesh.faces.push_back({base, base+1, base+2});
        }
        return mesh;
    }

    std::istringstream is(buf);
    std::string line, tok;
    int vc = 0;
    while (std::getline(is, line)) {
        std::istringstream ls(line);
        if (!(ls >> tok) || tok != "vertex") continue;
        double x, y, z;
        if (ls >> x >> y >> z) {
            mesh.vertices.push_back({x, y, z});
            if (++vc == 3) {
                const int b = (int)mesh.vertices.size() - 3;
                mesh.faces.push_back({b, b+1, b+2});
                vc = 0;
            }
        }
    }

    // Nothing parsed as ASCII, but the file carries a plausible binary header:
    // it is a truncated or padded binary STL. Say so, rather than letting the
    // caller see the misleading "contains no vertices".
    if (mesh.vertices.empty() && has_header && ntri > 0) {
        const unsigned long long want = 84ull + 50ull * ntri;
        throw std::runtime_error(
            "mesh '" + filename + "': truncated or padded binary STL - the header claims " +
            std::to_string(ntri) + " triangle(s), which needs " + std::to_string(want) +
            " bytes, but the file is " + std::to_string((unsigned long long)buf.size()) + " bytes");
    }
    return mesh;
}

bool ends_with_ci(const std::string& s, const char* ext) {
    const size_t n = std::strlen(ext);
    if (s.size() < n) return false;
    for (size_t i = 0; i < n; ++i)
        if (std::tolower((unsigned char)s[s.size()-n+i]) != ext[i]) return false;
    return true;
}

} // namespace

template <typename T>
Mesh<T> load_obj(const std::string& filename) {
    return build_mesh<T>(load_obj_raw(filename), filename);
}

template <typename T>
Mesh<T> load_stl(const std::string& filename) {
    return build_mesh<T>(load_stl_raw(filename), filename);
}

// Dispatch by file extension (.stl -> STL, otherwise OBJ).
template <typename T>
Mesh<T> load_mesh(const std::string& filename) {
    return ends_with_ci(filename, ".stl") ? load_stl<T>(filename) : load_obj<T>(filename);
}

template Mesh<float>  load_obj<float>(const std::string&);
template Mesh<double> load_obj<double>(const std::string&);
template Mesh<float>  load_stl<float>(const std::string&);
template Mesh<double> load_stl<double>(const std::string&);
template Mesh<float>  load_mesh<float>(const std::string&);
template Mesh<double> load_mesh<double>(const std::string&);

} // namespace flatland
