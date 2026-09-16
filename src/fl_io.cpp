#include "fl_io.hpp"

#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace flatland {

/* ----------------------
   IO & Loaders
   ---------------------- */
// Templated OBJ loader (vertices + triangulated faces) to support float/double.
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
    return mesh;
}

// STL loader, auto-detecting binary vs ASCII. STL stores independent per-triangle
// vertices (no shared indexing), so each triangle contributes three fresh vertices.
template <typename T>
Mesh<T> load_stl(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("cannot open mesh file '" + filename + "'");
    std::string buf((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    Mesh<T> mesh;

    // A binary STL is exactly 84 + 50*ntri bytes; this is the robust discriminator
    // (the leading "solid" keyword is NOT reliable — some binary files contain it).
    bool binary = false;
    uint32_t ntri = 0;
    if (buf.size() >= 84) {
        std::memcpy(&ntri, buf.data() + 80, 4);
        if (buf.size() == 84ull + 50ull * ntri) binary = true;
    }

    if (binary) {
        mesh.vertices.reserve((size_t)ntri * 3);
        mesh.faces.reserve(ntri);
        for (uint32_t t = 0; t < ntri; ++t) {
            const char* tri = buf.data() + 84 + (size_t)t * 50;
            float v[9];
            std::memcpy(v, tri + 12, 36);   // skip the 12-byte facet normal, read 3 vertices
            int base = (int)mesh.vertices.size();
            mesh.vertices.push_back({(T)v[0],(T)v[1],(T)v[2]});
            mesh.vertices.push_back({(T)v[3],(T)v[4],(T)v[5]});
            mesh.vertices.push_back({(T)v[6],(T)v[7],(T)v[8]});
            mesh.faces.push_back({base, base+1, base+2});
        }
    } else {
        std::istringstream is(buf);
        std::string line, tok;
        int vc = 0;
        while (std::getline(is, line)) {
            std::istringstream ls(line);
            if (!(ls >> tok) || tok != "vertex") continue;
            double x,y,z;
            if (ls >> x >> y >> z) {
                mesh.vertices.push_back({(T)x,(T)y,(T)z});
                if (++vc == 3) { int b=(int)mesh.vertices.size()-3; mesh.faces.push_back({b,b+1,b+2}); vc=0; }
            }
        }
    }
    return mesh;
}

// Dispatch by file extension (.stl -> STL, otherwise OBJ), then validate.
template <typename T>
Mesh<T> load_mesh(const std::string& filename) {
    auto ends_with_ci = [&](const char* ext) {
        size_t n = std::strlen(ext);
        if (filename.size() < n) return false;
        for (size_t i=0;i<n;++i)
            if (std::tolower(filename[filename.size()-n+i]) != ext[i]) return false;
        return true;
    };
    Mesh<T> mesh = ends_with_ci(".stl") ? load_stl<T>(filename) : load_obj<T>(filename);

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

template Mesh<float>  load_obj<float>(const std::string&);
template Mesh<double> load_obj<double>(const std::string&);
template Mesh<float>  load_stl<float>(const std::string&);
template Mesh<double> load_stl<double>(const std::string&);
template Mesh<float>  load_mesh<float>(const std::string&);
template Mesh<double> load_mesh<double>(const std::string&);

} // namespace flatland
