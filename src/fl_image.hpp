#pragma once

#include <string>
#include <cstdint>

#include "fl_raster.hpp"

namespace flatland {

// Write a false-colour PPM (P6) from raw buffers.
//
//   mask    w*h bytes, non-zero where a pixel is covered
//   vals    w*h values, or nullptr for a plain white silhouette
//   min_v / max_v   the range the colour ramp spans
//
// Row 0 of the buffers is the BOTTOM row in mesh space; PPM stores rows top
// down, so the writer emits them in reverse.
template <typename V>
void write_ppm(const std::string& fname, int w, int h,
               const uint8_t* mask, const V* vals, double min_v, double max_v);

// Convenience wrapper for a rendered view held in a Renderer.
template <typename T>
void save_ppm(const Renderer<T>& r, const std::string& fname, T min_v, T max_v, bool use_vals);

} // namespace flatland
