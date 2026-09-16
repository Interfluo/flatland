#pragma once

#include <string>
#include <vector>

#include "fl_vec.hpp"

namespace flatland {

struct Triangle3D {
    int v0_idx, v1_idx, v2_idx; // Indices into the vertex array
};

enum ValueMode { MODE_NONE, MODE_NODE, MODE_FACE };

// Read-only geometry, loaded once and shared across all views/threads.
//
// Vertices are stored RECENTERED: the bounding-box center of the input geometry
// is subtracted at load time and recorded in `origin`. Projection, rasterization
// and area are all translation-invariant in exact arithmetic, but not in float —
// a mesh sitting at a CAD or geospatial offset of 5e6 loses so much mantissa
// that its projected area was coming out ~25% high. Working in mesh-local
// coordinates removes that dependence on where the object happens to sit.
//
// `origin` exists so callers can be given back the frame they supplied;
// nothing in the pipeline itself needs it.
template <typename T>
struct Mesh {
    std::vector<Vec3<T>> vertices;
    std::vector<Triangle3D> faces;
    Vec3<double> origin{0,0,0};
};

// Geometry exactly as supplied — parsed from a file or handed over by a caller
// through the C API — before validation, recentering and narrowing. Always
// double: recentering has to happen before the conversion to the working type,
// or a distant mesh has already lost the mantissa bits that carry its detail.
struct RawMesh {
    std::vector<Vec3<double>> vertices;
    std::vector<Triangle3D> faces;
};

// Validate `raw`, recenter it on its bounding box, and narrow to T. `label`
// identifies the source in any error message. Throws std::runtime_error on
// empty geometry, an out-of-range index, or a non-finite coordinate.
template <typename T>
Mesh<T> build_mesh(RawMesh&& raw, const std::string& label);

// Narrow an already-built mesh, preserving its recentering.
template <typename To, typename From>
Mesh<To> narrow_mesh(const Mesh<From>& src) {
    Mesh<To> out;
    out.origin = src.origin;
    out.faces = src.faces;
    out.vertices.resize(src.vertices.size());
    for (size_t i = 0; i < src.vertices.size(); ++i)
        out.vertices[i] = { (To)src.vertices[i].x, (To)src.vertices[i].y, (To)src.vertices[i].z };
    return out;
}

} // namespace flatland
