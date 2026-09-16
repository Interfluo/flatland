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
template <typename T> inline Vec3<T> normalize(const Vec3<T>& v) { T n = std::sqrt(dot(v,v)); return n>0 ? Vec3<T>{v.x/n, v.y/n, v.z/n} : Vec3<T>{0,0,0}; }

// Convert an azimuth/elevation pair (degrees) into a unit view direction.
// Azimuth sweeps around +Z measured from +X; elevation rises above the XY plane.
inline Vec3<double> angle_to_dir(double az_deg, double el_deg) {
    const double PI = 3.14159265358979323846;
    double az = az_deg * PI / 180.0, el = el_deg * PI / 180.0;
    return { std::cos(el)*std::cos(az), std::cos(el)*std::sin(az), std::sin(el) };
}

} // namespace flatland
