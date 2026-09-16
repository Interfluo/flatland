#pragma once

#include <string>
#include <vector>

#include "fl_mesh.hpp"

namespace flatland {

struct BatchEntry {
    double nx, ny, nz;
    double resolution = -1.0; // -1 means use global default
    std::string data_file;    // empty means use global default
};

template <typename T>
void run_app(const std::string& obj, const std::string& default_data, const std::string& out_pre,
             const std::vector<BatchEntry>& batch, double default_res, bool cull, bool json,
             unsigned threads, ValueMode forced_mode);

} // namespace flatland
