#pragma once

#include <string>

#include "fl_field.hpp"
#include "fl_mesh.hpp"
#include "fl_raster.hpp"
#include "fl_vec.hpp"

namespace flatland {

template <typename T>
struct ViewResult {
    bool has_field = false;
    // True only when a field was supplied AND at least one pixel was covered.
    // When false the four statistics below are not measurements and must not be
    // reported as such; callers emit null rather than a fabricated zero.
    bool has_stats = false;
    T area = 0;
    T average_value = 0;
    T integral = 0;     // Σ value · pixel_area  (the area integral of the field)
    T min_val = 0;
    T max_val = 0;
    double time_seconds = 0.0;
    long covered_pixels = 0;
    int image_width = 0;
    int image_height = 0;
    std::string output_image;
};

template <typename T>
ViewResult<T> process_view(const Vec3<T>& view_dir, const Mesh<T>& mesh, const Field<T>* field,
                           T pix_sz, bool cull, Renderer<T>& renderer);

} // namespace flatland
