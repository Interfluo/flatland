// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Interfluo
//
// FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
// for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

#pragma once

#include <string>
#include <vector>
#include <cstddef>

#include "fl_mesh.hpp"

namespace flatland {

// A scalar field defined on the mesh. Kept separate from Mesh so different
// views can use different fields without mutating shared state (thread-safe).
template <typename T>
struct Field {
    std::vector<T> values;
    ValueMode mode = MODE_NONE;
};

/* ----------------------
   Scalar Fields
   ----------------------
   A field file is a matrix: one row per mesh entity (vertex or face), one or more
   whitespace-separated columns. A column is a timestep, so a whole time series can
   live in ONE file. The classic "one value per line" file is just the single-column
   case. A data token may select a column with a trailing '@<col>' (default 0). */

struct FieldToken { std::string path; size_t col; };

FieldToken parse_field_token(const std::string& token);

template <typename T>
struct FieldMatrix {
    std::vector<T> data;        // row-major: data[row*ncols + col]
    size_t nrows = 0, ncols = 0;
    ValueMode mode = MODE_NONE;
};

template <typename T>
void load_matrix_into(const std::string& filename, const Mesh<T>& mesh,
                      FieldMatrix<T>& m, ValueMode forced);

template <typename T>
void extract_column(const FieldMatrix<T>& m, size_t col, Field<T>& out);

} // namespace flatland
