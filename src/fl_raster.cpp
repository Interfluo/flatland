#include "fl_raster.hpp"

#include <algorithm>
#include <limits>
#include <cstddef>

namespace flatland {

template <typename T>
void Renderer<T>::resize(int w, int h, bool use_vals) {
    size_t needed = (size_t)w * h;
    if (needed == 0) { Nx = Ny = 0; return; }
    Nx = w; Ny = h;
    if (zbuffer.size() < needed) zbuffer.resize(needed);
    if (mask.size() < needed) mask.resize(needed);
    if (use_vals && val_buffer.size() < needed) val_buffer.resize(needed);

    std::fill(zbuffer.begin(), zbuffer.begin() + needed, std::numeric_limits<T>::max());
    std::fill(mask.begin(), mask.begin() + needed, 0);
    if (use_vals) std::fill(val_buffer.begin(), val_buffer.begin() + needed, 0);
}

template <typename T>
T edge_eval(const Vec2<T>& a, const Vec2<T>& b, const Vec2<T>& p) {
    return (p.x - a.x)*(b.y - a.y) - (p.y - a.y)*(b.x - a.x);
}

template <typename T>
void rasterize(const TriProjected<T>& tri, Renderer<T>& r, T xmin, T ymin, T pix_sz, bool use_vals) {
    T min_x = std::min({tri.v0.x, tri.v1.x, tri.v2.x});
    T max_x = std::max({tri.v0.x, tri.v1.x, tri.v2.x});
    T min_y = std::min({tri.v0.y, tri.v1.y, tri.v2.y});
    T max_y = std::max({tri.v0.y, tri.v1.y, tri.v2.y});

    int ix_min = std::max(0, static_cast<int>((min_x - xmin)/pix_sz));
    int ix_max = std::min(r.Nx - 1, static_cast<int>((max_x - xmin)/pix_sz));
    int iy_min = std::max(0, static_cast<int>((min_y - ymin)/pix_sz));
    int iy_max = std::min(r.Ny - 1, static_cast<int>((max_y - ymin)/pix_sz));

    Vec2<T> p_start = { xmin + (ix_min + (T)0.5)*pix_sz, ymin + (iy_min + (T)0.5)*pix_sz };

    // Edge function values at the first sampled pixel center
    T w0_row = edge_eval(tri.v1, tri.v2, p_start);
    T w1_row = edge_eval(tri.v2, tri.v0, p_start);
    T w2_row = edge_eval(tri.v0, tri.v1, p_start);

    // Per-pixel increments
    T A0 = (tri.v2.y - tri.v1.y)*pix_sz, B0 = (tri.v1.x - tri.v2.x)*pix_sz;
    T A1 = (tri.v0.y - tri.v2.y)*pix_sz, B1 = (tri.v2.x - tri.v0.x)*pix_sz;
    T A2 = (tri.v1.y - tri.v0.y)*pix_sz, B2 = (tri.v0.x - tri.v1.x)*pix_sz;

    for (int iy = iy_min; iy <= iy_max; ++iy) {
        T w0=w0_row, w1=w1_row, w2=w2_row;
        size_t idx = (size_t)iy * r.Nx + ix_min;
        for (int ix = ix_min; ix <= ix_max; ++ix) {
            if (w0 >= 0 && w1 >= 0 && w2 >= 0) {
                T z = (w0*tri.inv_area)*tri.z0 + (w1*tri.inv_area)*tri.z1 + (w2*tri.inv_area)*tri.z2;
                if (z < r.zbuffer[idx]) {
                    r.zbuffer[idx] = z;
                    r.mask[idx] = 1;
                    if (use_vals)
                        r.val_buffer[idx] = (w0*tri.val0 + w1*tri.val1 + w2*tri.val2)*tri.inv_area;
                }
            }
            w0+=A0; w1+=A1; w2+=A2; idx++;
        }
        w0_row+=B0; w1_row+=B1; w2_row+=B2;
    }
}

template struct TriProjected<float>;
template struct TriProjected<double>;
template class Renderer<float>;
template class Renderer<double>;

template float  edge_eval<float>(const Vec2<float>&, const Vec2<float>&, const Vec2<float>&);
template double edge_eval<double>(const Vec2<double>&, const Vec2<double>&, const Vec2<double>&);
template void rasterize<float>(const TriProjected<float>&, Renderer<float>&, float, float, float, bool);
template void rasterize<double>(const TriProjected<double>&, Renderer<double>&, double, double, double, bool);

} // namespace flatland
