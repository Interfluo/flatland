#pragma once

#include <string>

#include "fl_mesh.hpp"

namespace flatland {

template <typename T>
Mesh<T> load_obj(const std::string& filename);

template <typename T>
Mesh<T> load_stl(const std::string& filename);

template <typename T>
Mesh<T> load_mesh(const std::string& filename);

} // namespace flatland
