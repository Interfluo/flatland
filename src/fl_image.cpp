// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Interfluo
//
// FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
// for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

#include "fl_image.hpp"

#include "fl_deflate.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <cstring>

namespace flatland {

namespace {

const uint8_t BACKGROUND[3] = {30, 30, 35};

// A five-stop ramp: blue -> cyan -> yellow -> orange -> red.
void ramp_colour(double t, uint8_t rgb[3]) {
    static const double stops[5][3] = {{0,0,1}, {0,1,1}, {1,1,0}, {1,.5,0}, {1,0,0}};
    t = std::max(0.0, std::min(1.0, t)) * 4.0;
    int i = (int)t;
    if (i >= 4) i = 3;
    const double f = t - i;
    for (int c = 0; c < 3; ++c)
        rgb[c] = (uint8_t)(((1.0 - f) * stops[i][c] + f * stops[i+1][c]) * 255.0 + 0.5);
}

void put_u32(std::vector<uint8_t>& v, uint32_t x) {
    v.push_back((uint8_t)(x >> 24)); v.push_back((uint8_t)(x >> 16));
    v.push_back((uint8_t)(x >> 8));  v.push_back((uint8_t)x);
}

void png_chunk(std::ofstream& f, const char tag[4], const std::vector<uint8_t>& data) {
    std::vector<uint8_t> head;
    put_u32(head, (uint32_t)data.size());
    f.write((const char*)head.data(), (std::streamsize)head.size());
    f.write(tag, 4);
    if (!data.empty()) f.write((const char*)data.data(), (std::streamsize)data.size());

    // The CRC covers the type and the data, but not the length.
    uint32_t crc = crc32_bytes((const uint8_t*)tag, 4);
    // Continue the same CRC over the payload. crc32_bytes seeds with its
    // argument, so feeding the previous result chains the two runs.
    if (!data.empty()) crc = crc32_bytes(data.data(), data.size(), crc);
    std::vector<uint8_t> tail;
    put_u32(tail, crc);
    f.write((const char*)tail.data(), 4);
}

// Sum of the filtered bytes read as signed, the standard heuristic for picking
// a filter per scanline. Lower is more compressible.
uint64_t filter_cost(const std::vector<uint8_t>& row) {
    uint64_t sum = 0;
    for (uint8_t b : row) sum += (b < 128) ? b : (uint8_t)(256 - b);
    return sum;
}

uint8_t paeth(int a, int b, int c) {
    const int p = a + b - c;
    const int pa = std::abs(p - a), pb = std::abs(p - b), pc = std::abs(p - c);
    if (pa <= pb && pa <= pc) return (uint8_t)a;
    if (pb <= pc) return (uint8_t)b;
    return (uint8_t)c;
}

/*
 * Filter one scanline five ways and keep the cheapest. PNG allows the filter to
 * be chosen per row, and doing so is worth a good deal on rendered output: a
 * heat map has smooth horizontal gradients, where Sub and Paeth both leave long
 * runs of near-zero bytes for the compressor to eat.
 */
void filter_row(const uint8_t* cur, const uint8_t* prev, size_t stride, int bpp,
                std::vector<uint8_t>& out) {
    std::vector<uint8_t> cand[5];
    for (int t = 0; t < 5; ++t) cand[t].resize(stride);

    for (size_t i = 0; i < stride; ++i) {
        const int a = (i >= (size_t)bpp) ? cur[i - bpp] : 0;         // left
        const int b = prev ? prev[i] : 0;                             // above
        const int c = (prev && i >= (size_t)bpp) ? prev[i - bpp] : 0; // above-left
        const int x = cur[i];
        cand[0][i] = (uint8_t)x;
        cand[1][i] = (uint8_t)(x - a);
        cand[2][i] = (uint8_t)(x - b);
        cand[3][i] = (uint8_t)(x - ((a + b) >> 1));
        cand[4][i] = (uint8_t)(x - paeth(a, b, c));
    }

    int best = 0;
    uint64_t best_cost = std::numeric_limits<uint64_t>::max();
    for (int t = 0; t < 5; ++t) {
        const uint64_t cost = filter_cost(cand[t]);
        if (cost < best_cost) { best_cost = cost; best = t; }
    }
    out.push_back((uint8_t)best);
    out.insert(out.end(), cand[best].begin(), cand[best].end());
}

void check_raster(const std::string& fname, int w, int h, const uint8_t* mask) {
    if (w <= 0 || h <= 0)
        throw std::runtime_error("cannot write '" + fname + "': the view covered no pixels");
    if (!mask)
        throw std::runtime_error("cannot write '" + fname + "': no coverage mask");
}

} // namespace

/* ------------------------------------------------------------------------ */
/* PNG                                                                       */
/* ------------------------------------------------------------------------ */

template <typename V>
void write_png(const std::string& fname, int w, int h,
               const uint8_t* mask, const V* vals, double min_v, double max_v) {
    check_raster(fname, w, h, mask);

    // Normalize the ramp over the field's actual span. The flat-field guard is
    // RELATIVE: an absolute floor renders any field whose values are small,
    // however wide its dynamic range, as a single flat colour.
    double range = max_v - min_v;
    const double mag = std::max(std::abs(min_v), std::abs(max_v));
    if (!std::isfinite(range) || range <= 1e-12 * std::max(mag, 1e-300)) range = 1.0;

    const int bpp = 3;
    const size_t stride = (size_t)w * bpp;

    std::vector<uint8_t> raw;                  // filtered scanlines, top-down
    raw.reserve((stride + 1) * (size_t)h);
    std::vector<uint8_t> cur(stride), prev(stride);
    bool have_prev = false;

    for (int y = h - 1; y >= 0; --y) {         // flip: our row 0 is the bottom
        for (int x = 0; x < w; ++x) {
            const size_t i = (size_t)y * w + x;
            uint8_t* px = &cur[(size_t)x * bpp];
            if (!mask[i]) {
                std::memcpy(px, BACKGROUND, 3);
            } else if (!vals) {
                px[0] = px[1] = px[2] = 255;   // silhouette
            } else {
                ramp_colour(((double)vals[i] - min_v) / range, px);
            }
        }
        filter_row(cur.data(), have_prev ? prev.data() : nullptr, stride, bpp, raw);
        std::swap(cur, prev);
        have_prev = true;
    }

    const std::vector<uint8_t> idat = zlib_compress(raw.data(), raw.size(), 6);

    std::ofstream f(fname, std::ios::binary);
    if (!f.is_open()) throw std::runtime_error("cannot write image '" + fname + "'");

    static const uint8_t SIG[8] = {0x89, 'P', 'N', 'G', 0x0D, 0x0A, 0x1A, 0x0A};
    f.write((const char*)SIG, 8);

    std::vector<uint8_t> ihdr;
    put_u32(ihdr, (uint32_t)w);
    put_u32(ihdr, (uint32_t)h);
    ihdr.push_back(8);      // bit depth
    ihdr.push_back(2);      // colour type 2 = truecolour RGB
    ihdr.push_back(0);      // compression: deflate
    ihdr.push_back(0);      // filter method 0
    ihdr.push_back(0);      // no interlace
    png_chunk(f, "IHDR", ihdr);
    png_chunk(f, "IDAT", idat);
    png_chunk(f, "IEND", {});

    if (!f) throw std::runtime_error("failed while writing image '" + fname + "'");
}

/* ------------------------------------------------------------------------ */
/* NPY                                                                       */
/* ------------------------------------------------------------------------ */

template <typename V>
void write_npy(const std::string& fname, int w, int h,
               const uint8_t* mask, const V* vals) {
    check_raster(fname, w, h, mask);
    if (!vals)
        throw std::runtime_error("cannot write '" + fname +
                                 "': this view carried no field, so there are no values to export");

    // .npy v1.0: magic, version, 2-byte little-endian header length, then an
    // ASCII dict padded so the data starts on a 64-byte boundary.
    std::ostringstream dict;
    dict << "{'descr': '<f8', 'fortran_order': False, 'shape': (" << h << ", " << w << "), }";
    std::string header = dict.str();
    const size_t prelude = 10;                  // magic(6) + version(2) + len(2)
    size_t total = prelude + header.size() + 1; // +1 for the trailing newline
    const size_t padded = (total + 63) / 64 * 64;
    header.append(padded - total, ' ');
    header.push_back('\n');

    std::ofstream f(fname, std::ios::binary);
    if (!f.is_open()) throw std::runtime_error("cannot write array '" + fname + "'");

    f.write("\x93NUMPY", 6);
    const uint8_t ver[2] = {1, 0};
    f.write((const char*)ver, 2);
    const uint16_t hlen = (uint16_t)header.size();
    const uint8_t hl[2] = {(uint8_t)(hlen & 0xFF), (uint8_t)(hlen >> 8)};
    f.write((const char*)hl, 2);
    f.write(header.data(), (std::streamsize)header.size());

    // Row 0 stays the bottom row, matching every other FlatLand interface.
    // Uncovered pixels are NaN: nothing was measured there.
    const double nan_v = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> row((size_t)w);
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            const size_t i = (size_t)y * w + x;
            row[(size_t)x] = mask[i] ? (double)vals[i] : nan_v;
        }
        f.write((const char*)row.data(), (std::streamsize)(row.size() * sizeof(double)));
    }
    if (!f) throw std::runtime_error("failed while writing array '" + fname + "'");
}

/* ------------------------------------------------------------------------ */
/* Renderer wrappers                                                         */
/* ------------------------------------------------------------------------ */

template <typename T>
void save_png(const Renderer<T>& r, const std::string& fname, T min_v, T max_v, bool use_vals) {
    if (use_vals && r.val_buffer.size() < (size_t)r.Nx * r.Ny)
        throw std::runtime_error("cannot write image '" + fname + "': value buffer is smaller than the raster");
    write_png<T>(fname, r.Nx, r.Ny, r.mask.data(),
                 use_vals ? r.val_buffer.data() : nullptr,
                 (double)min_v, (double)max_v);
}

template <typename T>
void save_npy(const Renderer<T>& r, const std::string& fname, bool use_vals) {
    if (use_vals && r.val_buffer.size() < (size_t)r.Nx * r.Ny)
        throw std::runtime_error("cannot write array '" + fname + "': value buffer is smaller than the raster");
    write_npy<T>(fname, r.Nx, r.Ny, r.mask.data(),
                 use_vals ? r.val_buffer.data() : nullptr);
}

template void write_png<float>(const std::string&, int, int, const uint8_t*, const float*, double, double);
template void write_png<double>(const std::string&, int, int, const uint8_t*, const double*, double, double);
template void write_npy<float>(const std::string&, int, int, const uint8_t*, const float*);
template void write_npy<double>(const std::string&, int, int, const uint8_t*, const double*);

template void save_png<float>(const Renderer<float>&, const std::string&, float, float, bool);
template void save_png<double>(const Renderer<double>&, const std::string&, double, double, bool);
template void save_npy<float>(const Renderer<float>&, const std::string&, bool);
template void save_npy<double>(const Renderer<double>&, const std::string&, bool);

} // namespace flatland
