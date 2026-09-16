#pragma once

#include <vector>

#include "fl_vec.hpp"

namespace flatland {

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

} // namespace flatland
