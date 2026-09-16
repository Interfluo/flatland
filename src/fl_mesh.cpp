#include "fl_mesh.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>

namespace flatland {

template <typename T>
Mesh<T> build_mesh(RawMesh&& raw, const std::string& label) {
    if (raw.vertices.empty()) throw std::runtime_error("mesh '" + label + "' contains no vertices");
    if (raw.faces.empty())    throw std::runtime_error("mesh '" + label + "' contains no faces");

    const long long nv = (long long)raw.vertices.size();
    for (const auto& f : raw.faces) {
        if (f.v0_idx < 0 || f.v0_idx >= nv ||
            f.v1_idx < 0 || f.v1_idx >= nv ||
            f.v2_idx < 0 || f.v2_idx >= nv)
            throw std::runtime_error("mesh '" + label + "' references an out-of-range vertex index");
    }

    double lo[3] = { raw.vertices[0].x, raw.vertices[0].y, raw.vertices[0].z };
    double hi[3] = { lo[0], lo[1], lo[2] };
    for (const auto& p : raw.vertices) {
        if (!is_finite(p))
            throw std::runtime_error("mesh '" + label + "' contains a non-finite vertex coordinate");
        const double c[3] = { p.x, p.y, p.z };
        for (int k = 0; k < 3; ++k) { lo[k] = std::min(lo[k], c[k]); hi[k] = std::max(hi[k], c[k]); }
    }

    Mesh<T> mesh;
    // Recenter on the bounding-box center. Projection and area are
    // translation-invariant in exact arithmetic but not in float, so working in
    // mesh-local coordinates is what keeps a mesh at a CAD or geospatial offset
    // from reporting a materially wrong area.
    mesh.origin = { 0.5*(lo[0]+hi[0]), 0.5*(lo[1]+hi[1]), 0.5*(lo[2]+hi[2]) };
    if (!is_finite(mesh.origin))
        throw std::runtime_error("mesh '" + label + "' has a degenerate bounding box");

    mesh.vertices.resize(raw.vertices.size());
    for (size_t i = 0; i < raw.vertices.size(); ++i)
        mesh.vertices[i] = { (T)(raw.vertices[i].x - mesh.origin.x),
                             (T)(raw.vertices[i].y - mesh.origin.y),
                             (T)(raw.vertices[i].z - mesh.origin.z) };
    mesh.faces = std::move(raw.faces);
    return mesh;
}

template struct Mesh<float>;
template struct Mesh<double>;

template Mesh<float>  build_mesh<float>(RawMesh&&, const std::string&);
template Mesh<double> build_mesh<double>(RawMesh&&, const std::string&);

} // namespace flatland
