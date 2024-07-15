#ifndef POINT_H
#define POINT_H

#include <cmath>

#include "Types.h"

namespace HSDetection
{
    class Point
    {
    public:
        FloatGeom x;
        FloatGeom y;

        Point() : x(0), y(0) {}
        Point(FloatGeom x, FloatGeom y) : x(x), y(y) {}
        ~Point() {}

        FloatGeom abs() const { return sqrt(x * x + y * y); }

        Point &operator+=(const Point &other) { return x += other.x, y += other.y, *this; }
        Point &operator-=(const Point &other) { return x -= other.x, y -= other.y, *this; }
        Point &operator*=(FloatGeom mult) { return x *= mult, y *= mult, *this; }
        Point &operator/=(FloatGeom divisor) { return x /= divisor, y /= divisor, *this; }

        friend Point operator+(Point lhs, const Point &rhs) { return lhs += rhs; }
        friend Point operator-(Point lhs, const Point &rhs) { return lhs -= rhs; }
        friend Point operator*(Point pt, FloatGeom mult) { return pt *= mult; }
        friend Point operator*(FloatGeom mult, Point pt) { return pt *= mult; }
        friend Point operator/(Point pt, FloatGeom divisor) { return pt /= divisor; }
        // pass by value to allow optimization on chained ops
    };

} // namespace HSDetection

#endif
