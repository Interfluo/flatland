#include "fl_image.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <cstddef>
#include <cstdint>

namespace flatland {

/* ----------------------
   Image Output
   ---------------------- */
template <typename V>
void write_ppm(const std::string& fname, int w, int h,
               const uint8_t* mask, const V* vals, double min_v, double max_v) {
    if (w <= 0 || h <= 0)
        throw std::runtime_error("cannot write image '" + fname + "': the view covered no pixels");
    if (!mask)
        throw std::runtime_error("cannot write image '" + fname + "': no coverage mask");

    std::ofstream f(fname, std::ios::binary);
    if (!f.is_open()) throw std::runtime_error("cannot write image '" + fname + "'");
    f << "P6\n" << w << " " << h << "\n255\n";

    // Normalize the ramp over the field's actual span. The flat-field guard has
    // to be RELATIVE: an absolute floor renders any field whose values are
    // small — however wide its dynamic range — as a single flat colour.
    double range = max_v - min_v;
    const double mag = std::max(std::abs(min_v), std::abs(max_v));
    if (!std::isfinite(range) || range <= 1e-12 * std::max(mag, 1e-300)) range = 1.0;

    const uint8_t bg[3] = {30,30,35};
    // Simple 5-stop heatmap (blue -> cyan -> yellow -> orange -> red)
    const double stops[5][3] = {{0,0,1}, {0,1,1}, {1,1,0}, {1,.5,0}, {1,0,0}};

    for (int iy = h-1; iy >= 0; --iy) {
        for (int ix = 0; ix < w; ++ix) {
            const size_t i = (size_t)iy * w + ix;
            if (!mask[i]) { f.write((const char*)bg, 3); continue; }
            uint8_t rgb[3];
            if (!vals) {
                rgb[0]=rgb[1]=rgb[2]=255;            // silhouette: solid white
            } else {
                double t = ((double)vals[i] - min_v)/range;
                t = std::max(0.0, std::min(1.0, t)) * 4.0;
                int idx = (int)t;
                if (idx >= 4) idx = 3;
                const double fr = t - idx;
                for (int c=0;c<3;++c)
                    rgb[c] = (uint8_t)(((1.0-fr)*stops[idx][c] + fr*stops[idx+1][c])*255);
            }
            f.write((const char*)rgb, 3);
        }
    }
    if (!f) throw std::runtime_error("failed while writing image '" + fname + "'");
}

template <typename T>
void save_ppm(const Renderer<T>& r, const std::string& fname, T min_v, T max_v, bool use_vals) {
    if (use_vals && r.val_buffer.size() < (size_t)r.Nx * r.Ny)
        throw std::runtime_error("cannot write image '" + fname + "': value buffer is smaller than the raster");
    write_ppm<T>(fname, r.Nx, r.Ny, r.mask.data(),
                 use_vals ? r.val_buffer.data() : nullptr,
                 (double)min_v, (double)max_v);
}

template void write_ppm<float>(const std::string&, int, int, const uint8_t*, const float*, double, double);
template void write_ppm<double>(const std::string&, int, int, const uint8_t*, const double*, double, double);

template void save_ppm<float>(const Renderer<float>&, const std::string&, float, float, bool);
template void save_ppm<double>(const Renderer<double>&, const std::string&, double, double, bool);

} // namespace flatland
