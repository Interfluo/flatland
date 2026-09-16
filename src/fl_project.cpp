#include "fl_project.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <vector>
#include <cstddef>

namespace flatland {

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

    // Orthonormal view basis (n = view normal, u/v = screen axes)
    Vec3<T> n = normalize(view_dir);
    Vec3<T> h = (std::abs(n.x) < 0.9) ? Vec3<T>{1,0,0} : Vec3<T>{0,1,0};
    Vec3<T> u = normalize(cross(h, n));
    Vec3<T> v = cross(n, u);

    // Project triangles to the view plane
    std::vector<TriProjected<T>> tris;
    tris.reserve(mesh.faces.size());
    T b_min_x=1e20, b_max_x=-1e20, b_min_y=1e20, b_max_y=-1e20;

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

        T area2 = edge_eval(p0, p1, p2); // twice signed screen-space area; >0 = CCW
        if (std::abs(area2) < 1e-12) continue;

        b_min_x = std::min({b_min_x, p0.x, p1.x, p2.x});
        b_max_x = std::max({b_max_x, p0.x, p1.x, p2.x});
        b_min_y = std::min({b_min_y, p0.y, p1.y, p2.y});
        b_max_y = std::max({b_max_y, p0.y, p1.y, p2.y});

        // Normalize winding to CCW so edge functions are positive inside.
        if (area2 < 0) {
            tris.push_back({ p0, p2, p1,
                             dot(v0,n), dot(v2,n), dot(v1,n),
                             val0, val2, val1,
                             (T)(1.0/(-area2)) });
        } else {
            tris.push_back({ p0, p1, p2,
                             dot(v0,n), dot(v1,n), dot(v2,n),
                             val0, val1, val2,
                             (T)(1.0/area2) });
        }
    }

    if (!tris.empty()) {
        b_min_x -= pix_sz; b_max_x += pix_sz;
        b_min_y -= pix_sz; b_max_y += pix_sz;
        int w = (int)std::ceil((b_max_x - b_min_x)/pix_sz);
        int h = (int)std::ceil((b_max_y - b_min_y)/pix_sz);

        renderer.resize(w, h, use_vals);
        res.image_width=w; res.image_height=h;

        for (const auto& t : tris) rasterize(t, renderer, b_min_x, b_min_y, pix_sz, use_vals);

        // Aggregate statistics over visible pixels
        size_t npix = (size_t)w * h;
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
        T pix_area = pix_sz * pix_sz;
        res.area = res.covered_pixels * pix_area;
        if (use_vals && res.covered_pixels) {
            res.average_value = (T)(sum / res.covered_pixels);
            res.integral = (T)(sum * pix_area);  // ∫ f dA over the visible projection
            res.min_val = mn;
            res.max_val = mx;
        }
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
