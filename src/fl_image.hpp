// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Interfluo
//
// FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
// for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

#pragma once

#include <string>
#include <cstdint>

#include "fl_raster.hpp"

namespace flatland {

/*
 * Raster output.
 *
 * Two kinds, because they answer different questions:
 *
 *   PNG  a picture. The field is mapped through a colour ramp to 8 bits per
 *        channel, which is right for looking at and wrong for computing with.
 *   NPY  the numbers. The interpolated field value at every pixel, as float64,
 *        with NaN where nothing was covered. This is what you want if the
 *        raster is an input to something else rather than a figure.
 *
 * Both take the raster BOTTOM-UP, row 0 first, matching the rest of FlatLand.
 * PNG stores rows top-down, so write_png flips on the way out; NPY keeps the
 * bottom-up order, so a caller plotting it wants origin='lower'.
 */

// False-colour PNG. `vals` may be null, giving a white silhouette on the
// background colour. min_v/max_v set the ends of the ramp.
template <typename V>
void write_png(const std::string& fname, int w, int h,
               const uint8_t* mask, const V* vals, double min_v, double max_v);

// NumPy .npy, float64, shape (h, w), C order, bottom-up. Uncovered pixels are
// NaN rather than zero: nothing was measured there, and zero is a value a field
// can legitimately take.
template <typename V>
void write_npy(const std::string& fname, int w, int h,
               const uint8_t* mask, const V* vals);

// Convenience wrappers for a rendered view held in a Renderer.
template <typename T>
void save_png(const Renderer<T>& r, const std::string& fname, T min_v, T max_v, bool use_vals);

template <typename T>
void save_npy(const Renderer<T>& r, const std::string& fname, bool use_vals);

} // namespace flatland
