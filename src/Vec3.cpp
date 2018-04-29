#include <cassert>
#include <cmath>
#include "Vec3.hpp"

Vec3::Vec3()
{
}

Vec3::Vec3(const double x, const double y, const double z)
{
    n[0] = x;
    n[1] = y;
    n[2] = z;
}

Vec3::Vec3(const Vec3 &v)
{
    n[0] = v.n[0];
    n[1] = v.n[1];
    n[2] = v.n[2];
}

// ASSIGNMENT OPERATORS

Vec3 &Vec3::operator=(const Vec3 &v)
{
    n[0] = v.n[0];
    n[1] = v.n[1];
    n[2] = v.n[2];
    return *this;
}

Vec3 &Vec3::operator+=(const Vec3 &v)
{
    n[0] += v.n[0];
    n[1] += v.n[1];
    n[2] += v.n[2];
    return *this;
}

Vec3 &Vec3::operator-=(const Vec3 &v)
{
    n[0] -= v.n[0];
    n[1] -= v.n[1];
    n[2] -= v.n[2];
    return *this;
}

Vec3 &Vec3::operator*=(const double d)
{
    n[0] *= d;
    n[1] *= d;
    n[2] *= d;
    return *this;
}

Vec3 &Vec3::operator/=(const double d)
{
    double d_inv = 1.0f / d;
    n[0] *= d_inv;
    n[1] *= d_inv;
    n[2] *= d_inv;
    return *this;
}

double &Vec3::operator[](int i)
{
    assert(!(i < 0 || i > 2));
    return n[i];
}

double Vec3::operator[](int i) const
{
    assert(!(i < 0 || i > 2));
    return n[i];
}

void Vec3::set(const double x, const double y, const double z)
{
    n[0] = x;
    n[1] = y;
    n[2] = z;
}

// SPECIAL FUNCTIONS

double Vec3::norm() const
{
    return std::sqrt(sqrLength());
}

double Vec3::sqrLength() const
{
    return n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
}

Vec3 &Vec3::normalize() // it is up to caller to avoid divide-by-zero
{
    double len = norm();
    if (len > 0.000001)
        *this /= norm();
    return *this;
}

Vec3 Vec3::cross(const Vec3 &v) const
{
    Vec3 tmp;
    tmp[0] = n[1] * v.n[2] - n[2] * v.n[1];
    tmp[1] = n[2] * v.n[0] - n[0] * v.n[2];
    tmp[2] = n[0] * v.n[1] - n[1] * v.n[0];
    return tmp;
}

Vec3 operator-(const Vec3 &a)
{
    return Vec3(-a.n[0], -a.n[1], -a.n[2]);
}

Vec3 operator+(const Vec3 &a, const Vec3 &b)
{
    return Vec3(a.n[0] + b.n[0], a.n[1] + b.n[1], a.n[2] + b.n[2]);
}

Vec3 operator-(const Vec3 &a, const Vec3 &b)
{
    return Vec3(a.n[0] - b.n[0], a.n[1] - b.n[1], a.n[2] - b.n[2]);
}

Vec3 operator*(const Vec3 &a, const double d)
{
    return Vec3(d * a.n[0], d * a.n[1], d * a.n[2]);
}

Vec3 operator*(const double d, const Vec3 &a)
{
    return a * d;
}

Vec3 operator*(const Vec3 &a, const Vec3 &b)
{
    return Vec3(a.n[0] * b.n[0], a.n[1] * b.n[1], a.n[2] * b.n[2]);
}

Vec3 operator/(const Vec3 &a, const double d)
{
    double d_inv = 1.0f / d;
    return Vec3(a.n[0] * d_inv, a.n[1] * d_inv, a.n[2] * d_inv);
}

Vec3 operator^(const Vec3 &a, const Vec3 &b)
{
    return Vec3(a.n[1] * b.n[2] - a.n[2] * b.n[1],
                a.n[2] * b.n[0] - a.n[0] * b.n[2],
                a.n[0] * b.n[1] - a.n[1] * b.n[0]);
}

bool operator==(const Vec3 &a, const Vec3 &b)
{
    return (a.n[0] == b.n[0]) && (a.n[1] == b.n[1]) && (a.n[2] == b.n[2]);
}

bool operator!=(const Vec3 &a, const Vec3 &b)
{
    return !(a == b);
}

Vec3 prod(const Vec3 &a, const Vec3 &b)
{
    return Vec3(a.n[0] * b.n[0], a.n[1] * b.n[1], a.n[2] * b.n[2]);
}

double dot(const Vec3 &a, const Vec3 &b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double distance(const Vec3 &a, const Vec3 &b) // distance
{
    return std::sqrt((b[0] - a[0]) * (b[0] - a[0]) +
                     (b[1] - a[1]) * (b[1] - a[1]) +
                     (b[2] - a[2]) * (b[2] - a[2]));
}

double distanceSqr(const Vec3 &a, const Vec3 &b) // distance
{
    return ((b[0] - a[0]) * (b[0] - a[0]) +
            (b[1] - a[1]) * (b[1] - a[1]) +
            (b[2] - a[2]) * (b[2] - a[2]));
}