#pragma once

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

} // namespace flatland
