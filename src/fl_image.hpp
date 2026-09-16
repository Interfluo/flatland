#pragma once

#include <string>

#include "fl_raster.hpp"

namespace flatland {

template <typename T>
void save_ppm(const Renderer<T>& r, const std::string& fname, T min_v, T max_v, bool use_vals);

} // namespace flatland
