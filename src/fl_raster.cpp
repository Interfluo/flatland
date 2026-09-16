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

    const T px0 = xmin + (ix_min + (T)0.5)*pix_sz;

    // Per-pixel increments along a row.
    const T A0 = (tri.v2.y - tri.v1.y)*pix_sz;
    const T A1 = (tri.v0.y - tri.v2.y)*pix_sz;
    const T A2 = (tri.v1.y - tri.v0.y)*pix_sz;

    // Interpolate RELATIVE to vertex 0 rather than as a weighted sum of all
    // three vertices. The two are equivalent in exact arithmetic because
    // w0+w1+w2 == area2, but the edge functions are stepped incrementally while
    // inv_area is fixed, so in float that identity drifts and the weights stop
    // summing to one. Expressed this way the drift multiplies the DIFFERENCES
    // between vertex values, so a constant field interpolates to exactly that
    // constant and the error on a varying field scales with its range rather
    // than its magnitude. Previously a constant field of 1e7 came back spanning
    // 9999999 to 10000733.
    const T dz1 = tri.z1 - tri.z0, dz2 = tri.z2 - tri.z0;
    const T dv1 = tri.val1 - tri.val0, dv2 = tri.val2 - tri.val0;

    for (int iy = iy_min; iy <= iy_max; ++iy) {
        // Recompute each row's starting edge values exactly instead of carrying
        // an accumulator down the rows: the error in a running sum grows with
        // the raster height, which is precisely when accuracy matters most.
        const Vec2<T> p_row = { px0, ymin + (iy + (T)0.5)*pix_sz };
        T w0 = edge_eval(tri.v1, tri.v2, p_row);
        T w1 = edge_eval(tri.v2, tri.v0, p_row);
        T w2 = edge_eval(tri.v0, tri.v1, p_row);

        size_t idx = (size_t)iy * r.Nx + ix_min;
        for (int ix = ix_min; ix <= ix_max; ++ix) {
            if (w0 >= 0 && w1 >= 0 && w2 >= 0) {
                const T z = tri.z0 + (w1*dz1 + w2*dz2)*tri.inv_area;
                if (z < r.zbuffer[idx]) {
                    r.zbuffer[idx] = z;
                    r.mask[idx] = 1;
                    if (use_vals)
                        r.val_buffer[idx] = tri.val0 + (w1*dv1 + w2*dv2)*tri.inv_area;
                }
            }
            w0+=A0; w1+=A1; w2+=A2; idx++;
        }
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
