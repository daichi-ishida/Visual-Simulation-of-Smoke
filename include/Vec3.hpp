#pragma once

class Vec3
{
public:
  double n[3];

  // Constructors
  Vec3();
  Vec3(const double x, const double y, const double z);
  Vec3(const Vec3 &v); // copy constructor

  // Assignment operators
  Vec3 &operator=(const Vec3 &v);   // assignment of a Vec3
  Vec3 &operator+=(const Vec3 &v);  // incrementation by a Vec3
  Vec3 &operator-=(const Vec3 &v);  // decrementation by a Vec3
  Vec3 &operator*=(const double d); // multiplication by a constant
  Vec3 &operator/=(const double d); // division by a constant
  double &operator[](int i);        // indexing
  double operator[](int i) const;   // read-only indexing
  void set(const double x, const double y, const double z);

  // special functions
  double norm() const;             // length of a Vec3
  double sqrLength() const;        // squared length of a Vec3
  Vec3 &normalize();               // normalize a Vec3 in place
  Vec3 cross(const Vec3 &v) const; // cross product: self cross v

  // friends
  friend Vec3 operator-(const Vec3 &v);                    // -v1
  friend Vec3 operator+(const Vec3 &a, const Vec3 &b);     // v1 + v2
  friend Vec3 operator-(const Vec3 &a, const Vec3 &b);     // v1 - v2
  friend Vec3 operator*(const Vec3 &a, const double d);    // v1 * scalar
  friend Vec3 operator*(const double d, const Vec3 &a);    // scalar * v1
  friend Vec3 operator*(const Vec3 &a, const Vec3 &b);     // piecewise muliply
  friend Vec3 operator/(const Vec3 &a, const double d);    // v1 / scalar
  friend Vec3 operator^(const Vec3 &a, const Vec3 &b);     // cross product
  friend bool operator==(const Vec3 &a, const Vec3 &b);    // v1 == v2 ?
  friend bool operator!=(const Vec3 &a, const Vec3 &b);    // v1 != v2 ?
  friend Vec3 Prod(const Vec3 &a, const Vec3 &b);          // term by term *
  friend double Dot(const Vec3 &a, const Vec3 &b);         // dot product
  friend double Distance(const Vec3 &a, const Vec3 &b);    // distance
  friend double DistanceSqr(const Vec3 &a, const Vec3 &b); // distance sqr
};