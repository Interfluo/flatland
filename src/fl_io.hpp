// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Interfluo
//
// FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
// for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

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
