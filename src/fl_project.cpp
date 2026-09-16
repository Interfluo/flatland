#include "fl_project.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstddef>

namespace flatland {

namespace {

// Ceiling on the raster a single view may allocate. A legitimate view is far
// below this; passing a resolution finer than the geometry warrants is the only
// realistic way to exceed it, and the user deserves to be told that in those
// terms rather than through a signed-overflow crash or an allocator exception.
const double MAX_RASTER_PIXELS = 1.0e9;

} // namespace

/* ----------------------
   Processing Pipeline
   ---------------------- */
template <typename T>
ViewResult<T> process_view(const Vec3<T>& view_dir, const Mesh<T>& mesh, const Field<T>* field,
                           T pix_sz, bool cull, Renderer<T>& renderer) {
    auto t0 = std::chrono::high_resolution_clock::now();
    ViewResult<T> res;
    const bool use_vals = (field && field->mode != MODE_NONE);
    res.has_field = use_vals;

    if (!(pix_sz > 0) || !std::isfinite((double)pix_sz))
        throw std::runtime_error("resolution must be a positive, finite number");

    // Orthonormal view basis (n = view normal, u/v = screen axes). normalize()
    // returns the zero vector for a zero or non-finite direction.
    Vec3<T> n = normalize(view_dir);
    if (n.x == 0 && n.y == 0 && n.z == 0)
        throw std::runtime_error("view direction must be a non-zero, finite vector");
    Vec3<T> h = (std::abs(n.x) < 0.9) ? Vec3<T>{1,0,0} : Vec3<T>{0,1,0};
    Vec3<T> u = normalize(cross(h, n));
    Vec3<T> v = cross(n, u);

    // Project triangles to the view plane.
    std::vector<TriProjected<T>> tris;
    tris.reserve(mesh.faces.size());

    // Seed the bounds from the type, not a magic 1e20: a mesh whose projected
    // extent exceeded that constant used to end up with a bounding box clamped
    // to the sentinel, and a raster sized from the wrong span.
    bool have_bounds = false;
    T b_min_x = 0, b_max_x = 0, b_min_y = 0, b_max_y = 0;

    for (size_t i=0; i<mesh.faces.size(); ++i) {
        const auto& f = mesh.faces[i];
        const Vec3<T>& v0 = mesh.vertices[f.v0_idx];
        const Vec3<T>& v1 = mesh.vertices[f.v1_idx];
        const Vec3<T>& v2 = mesh.vertices[f.v2_idx];

        // The camera looks ALONG n, so it sits at -infinity on the n axis and the
        // z-buffer below keeps the smallest dot(v,n) as nearest. A front-facing
        // triangle therefore has an outward normal opposing n: dot(tri_n, n) < 0.
        // Culling must agree with that depth convention, otherwise --no-cull and the
        // default path resolve to opposite surfaces of a closed mesh.
        Vec3<T> tri_n = cross(v1-v0, v2-v0);
        if (cull && dot(tri_n, n) >= 0) continue;   // drop faces pointing away from the viewer

        Vec2<T> p0 = {dot(v0,u), dot(v0,v)};
        Vec2<T> p1 = {dot(v1,u), dot(v1,v)};
        Vec2<T> p2 = {dot(v2,u), dot(v2,v)};

        T val0=0, val1=0, val2=0;
        if (use_vals) {
            if (field->mode == MODE_NODE) {
                val0=field->values[f.v0_idx]; val1=field->values[f.v1_idx]; val2=field->values[f.v2_idx];
            } else { // MODE_FACE
                val0=val1=val2=field->values[i];
            }
        }

        // Twice the signed screen-space area. Note edge_eval's sign convention is
        // the opposite of the usual 2D cross product, so a CCW triangle gives a
        // NEGATIVE value here; the normalization below makes the stored winding
        // consistent either way.
        const T area2 = edge_eval(p0, p1, p2);
        if (!std::isfinite((double)area2)) continue;

        // Reject degeneracy RELATIVE to the triangle's own size. A fixed epsilon
        // is meaningless: it discards a legitimately small mesh outright, and at
        // ordinary coordinate magnitudes float roundoff dwarfs it, so genuinely
        // collinear triangles slip through and rasterize.
        const T ex1 = p1.x-p0.x, ey1 = p1.y-p0.y;
        const T ex2 = p2.x-p0.x, ey2 = p2.y-p0.y;
        const T ex3 = p2.x-p1.x, ey3 = p2.y-p1.y;
        const T scale = std::max({ ex1*ex1+ey1*ey1, ex2*ex2+ey2*ey2, ex3*ex3+ey3*ey3 });
        const T degenerate_tol = (T)(8 * (double)std::numeric_limits<T>::epsilon()) * scale;
        if (!(std::abs(area2) > degenerate_tol)) continue;

        if (!have_bounds) {
            b_min_x = b_max_x = p0.x;
            b_min_y = b_max_y = p0.y;
            have_bounds = true;
        }
        b_min_x = std::min({b_min_x, p0.x, p1.x, p2.x});
        b_max_x = std::max({b_max_x, p0.x, p1.x, p2.x});
        b_min_y = std::min({b_min_y, p0.y, p1.y, p2.y});
        b_max_y = std::max({b_max_y, p0.y, p1.y, p2.y});

        // Normalize winding so the edge functions are positive inside.
        if (area2 < 0) {
            tris.push_back({ p0, p2, p1,
                             dot(v0,n), dot(v2,n), dot(v1,n),
                             val0, val2, val1,
                             (T)(1.0/(-(double)area2)) });
        } else {
            tris.push_back({ p0, p1, p2,
                             dot(v0,n), dot(v1,n), dot(v2,n),
                             val0, val1, val2,
                             (T)(1.0/(double)area2) });
        }
    }

    if (tris.empty()) {
        // No geometry reached the raster: an entirely back-facing view, or one
        // seen exactly edge-on. The renderer is REUSED across views by its
        // worker, so it must be cleared rather than left holding the previous
        // view's dimensions and buffers.
        renderer.resize(0, 0, false);
        res.time_seconds = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0).count();
        return res;
    }

    b_min_x -= pix_sz; b_max_x += pix_sz;
    b_min_y -= pix_sz; b_max_y += pix_sz;

    // Size the raster in double and range-check before narrowing: (int)ceil of
    // an overlarge or non-finite quotient is undefined behavior, and in practice
    // lands on INT_MIN and then surfaces as an allocator error.
    const double wf = std::ceil(((double)b_max_x - (double)b_min_x)/(double)pix_sz);
    const double hf = std::ceil(((double)b_max_y - (double)b_min_y)/(double)pix_sz);
    if (!std::isfinite(wf) || !std::isfinite(hf) || wf < 1 || hf < 1 || wf*hf > MAX_RASTER_PIXELS) {
        std::ostringstream msg;
        msg << "resolution " << (double)pix_sz << " is too fine for this view: it needs a "
            << wf << " x " << hf << " raster (" << wf*hf << " pixels, limit "
            << MAX_RASTER_PIXELS << "). Use a larger -r / resolution.";
        throw std::runtime_error(msg.str());
    }
    const int rw = (int)wf, rh = (int)hf;

    try {
        renderer.resize(rw, rh, use_vals);
    } catch (const std::bad_alloc&) {
        throw std::runtime_error("out of memory allocating a " + std::to_string(rw) + " x " +
                                 std::to_string(rh) + " raster; use a larger resolution");
    }
    res.image_width=rw; res.image_height=rh;

    for (const auto& t : tris) rasterize(t, renderer, b_min_x, b_min_y, pix_sz, use_vals);

    // Aggregate statistics over visible pixels
    const size_t npix = (size_t)rw * rh;
    double sum = 0;
    T mn = std::numeric_limits<T>::max();
    T mx = -std::numeric_limits<T>::max();
    for (size_t k=0; k<npix; ++k) {
        if (!renderer.mask[k]) continue;
        res.covered_pixels++;
        if (use_vals) {
            T val = renderer.val_buffer[k];
            sum += val;
            if (val < mn) mn = val;
            if (val > mx) mx = val;
        }
    }
    const T pix_area = pix_sz * pix_sz;
    res.area = res.covered_pixels * pix_area;

    // The field statistics exist only if something was actually covered. When
    // nothing was, they are left at zero and has_stats says so: reporting a
    // fabricated 0 would corrupt any min/max taken across a batch.
    if (use_vals && res.covered_pixels > 0) {
        res.has_stats = true;
        res.average_value = (T)(sum / res.covered_pixels);
        res.integral = (T)(sum * pix_area);  // ∫ f dA over the visible projection
        res.min_val = mn;
        res.max_val = mx;
    }

    res.time_seconds = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0).count();
    return res;
}

template struct ViewResult<float>;
template struct ViewResult<double>;

template ViewResult<float>  process_view<float>(const Vec3<float>&, const Mesh<float>&, const Field<float>*,
                                                float, bool, Renderer<float>&);
template ViewResult<double> process_view<double>(const Vec3<double>&, const Mesh<double>&, const Field<double>*,
                                                 double, bool, Renderer<double>&);

} // namespace flatland
