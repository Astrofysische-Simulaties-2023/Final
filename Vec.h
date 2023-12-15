#ifndef Q
#define Q

#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

class Vec {
    double _x;
    double _y;
    double _z;
public:
    Vec() { _x=0; _y=0; _z=0;}
    Vec(double x, double y, double z)
        { _x = x; _y = y; _z = z;}
    double x() const { return _x; }
    double y() const { return _y; }
    double z() const { return _z; }

    double norm() const
        { return sqrt(_x*_x + _y*_y+ _z*_z); }
    double norm2() const
        { return _x*_x + _y*_y+ _z*_z;}
    double norm3() const
        { double r = sqrt(_x*_x + _y*_y + _z*_z); return r*r*r; }

    Vec& operator+=(Vec v) {
        _x += v._x;
        _y += v._y;
        _z += v._z;
        return *this;
    }

    Vec& operator-=(Vec v) {
        _x -= v._x;
        _y -= v._y;
        _z -= v._z;
        return *this;
    }
    Vec& operator*=(double s) {
        _x *= s;
        _y *= s;
        _z *= s;
        return *this;
    }
    Vec& operator/=(double s) {
        _x /= s;
        _y /= s;
        _z /= s;
        return *this;
    }

};

Vec operator+(Vec a, Vec b) { return a += b; }
Vec operator-(Vec a, Vec b) { return a -= b; }
Vec operator*(Vec a, double s) { return a *= s; }
Vec operator*(double s, Vec b) { return b *= s; }
Vec operator/(Vec a, double s) { return a /= s; }

#endif