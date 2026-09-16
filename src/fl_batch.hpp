// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Interfluo
//
// FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
// for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

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
             const std::string& npy_pre, const std::vector<BatchEntry>& batch,
             double default_res, bool cull, bool json, unsigned threads, ValueMode forced_mode);

} // namespace flatland
