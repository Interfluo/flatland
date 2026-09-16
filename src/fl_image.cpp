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
template <typename T>
void save_ppm(const Renderer<T>& r, const std::string& fname, T min_v, T max_v, bool use_vals) {
    std::ofstream f(fname, std::ios::binary);
    if (!f.is_open()) throw std::runtime_error("cannot write image '" + fname + "'");
    f << "P6\n" << r.Nx << " " << r.Ny << "\n255\n";
    double range = (double)(max_v - min_v);
    if (std::abs(range) < 1e-9) range = 1.0;
    const uint8_t bg[3] = {30,30,35};

    // Simple 5-stop heatmap (blue -> cyan -> yellow -> orange -> red)
    const double stops[5][3] = {{0,0,1}, {0,1,1}, {1,1,0}, {1,.5,0}, {1,0,0}};

    for (int iy=r.Ny-1; iy>=0; --iy) {
        for (int ix=0; ix<r.Nx; ++ix) {
            size_t i = (size_t)iy * r.Nx + ix;
            if (!r.mask[i]) { f.write((const char*)bg, 3); continue; }
            uint8_t rgb[3];
            if (!use_vals) {
                rgb[0]=rgb[1]=rgb[2]=255;            // silhouette: solid white
            } else {
                double t = ((double)r.val_buffer[i] - min_v)/range;
                t = std::max(0.0, std::min(1.0, t)) * 4.0;
                int idx = (int)t;
                if (idx >= 4) idx = 3;
                double fr = t - idx;
                for(int c=0;c<3;++c) rgb[c] = (uint8_t)(((1.0-fr)*stops[idx][c] + fr*stops[idx+1][c])*255);
            }
            f.write((const char*)rgb, 3);
        }
    }
}

template void save_ppm<float>(const Renderer<float>&, const std::string&, float, float, bool);
template void save_ppm<double>(const Renderer<double>&, const std::string&, double, double, bool);

} // namespace flatland
