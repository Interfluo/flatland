#pragma once

#include <cmath>

namespace flatland {

/* ----------------------
   Math & Structs (Templated)
   ---------------------- */
template <typename T> struct Vec3 { T x, y, z; };
template <typename T> struct Vec2 { T x, y; };

/* ----------------------
   Math Ops
   ---------------------- */
template <typename T> inline T dot(const Vec3<T>& a, const Vec3<T>& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
template <typename T> inline Vec3<T> cross(const Vec3<T>& a, const Vec3<T>& b) { return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x}; }
template <typename T> inline Vec3<T> operator-(const Vec3<T>& a, const Vec3<T>& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
// A direction's magnitude is irrelevant, so normalizing must not care about it
// either. Computing dot(v,v) in T lets an ordinary-looking direction such as
// (0,0,1e25) overflow to infinity in float — and (0,0,1e-25) flush to zero —
// after which the caller silently gets an empty projection. Scaling by the
// largest component first, and accumulating in double regardless of T, makes the
// result depend only on the direction. Returns the zero vector for a zero or
// non-finite input; callers treat that as an error.
template <typename T> inline Vec3<T> normalize(const Vec3<T>& v) {
    const double x = (double)v.x, y = (double)v.y, z = (double)v.z;
    const double m = std::fmax(std::fmax(std::fabs(x), std::fabs(y)), std::fabs(z));
    if (!(m > 0) || !std::isfinite(m)) return Vec3<T>{0,0,0};
    const double sx = x/m, sy = y/m, sz = z/m;
    const double n = std::sqrt(sx*sx + sy*sy + sz*sz);
    if (!(n > 0)) return Vec3<T>{0,0,0};
    return Vec3<T>{ (T)(sx/n), (T)(sy/n), (T)(sz/n) };
}

// True when every component is finite. A NaN or infinity reaching the rasterizer
// corrupts the bounding box and then the raster dimensions.
template <typename T> inline bool is_finite(const Vec3<T>& v) {
    return std::isfinite((double)v.x) && std::isfinite((double)v.y) && std::isfinite((double)v.z);
}

// Convert an azimuth/elevation pair (degrees) into a unit view direction.
// Azimuth sweeps around +Z measured from +X; elevation rises above the XY plane.
inline Vec3<double> angle_to_dir(double az_deg, double el_deg) {
    const double PI = 3.14159265358979323846;
    double az = az_deg * PI / 180.0, el = el_deg * PI / 180.0;
    return { std::cos(el)*std::cos(az), std::cos(el)*std::sin(az), std::sin(el) };
}

} // namespace flatland
