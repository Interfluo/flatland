#pragma once

#include <vector>
#include <cstdint>

#include "fl_vec.hpp"

namespace flatland {

template <typename T>
struct TriProjected {
    Vec2<T> v0, v1, v2;
    T z0, z1, z2;
    T val0, val1, val2;
    T inv_area;
};

/* ----------------------
   Renderer (one per worker thread)
   ---------------------- */
template <typename T>
class Renderer {
public:
    int Nx = 0, Ny = 0;
    std::vector<T> zbuffer;
    std::vector<T> val_buffer;
    std::vector<uint8_t> mask;

    void resize(int w, int h, bool use_vals);
};

/* ----------------------
   Rasterization Logic
   ---------------------- */
template <typename T>
T edge_eval(const Vec2<T>& a, const Vec2<T>& b, const Vec2<T>& p);

template <typename T>
void rasterize(const TriProjected<T>& tri, Renderer<T>& r, T xmin, T ymin, T pix_sz, bool use_vals);

} // namespace flatland
