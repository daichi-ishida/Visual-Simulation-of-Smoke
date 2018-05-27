#pragma once
#include <array>
#include "constants.hpp"
#include "Vec3.hpp"

template <int X, int Y, int Z>
class GridData
{
public:
  GridData();
  ~GridData();

  double &operator()(int i, int j, int k);
  double *begin();
  double *end();

  double interp(const Vec3 &pt);

private:
  double linearInterpolation(const Vec3 &pt);
  double monotonicCubicInterpolation(const Vec3 &pt);
  double axis_monotonicCubicInterpolation(const double f[], const double t) const;
  int sign(const double a) const;
  int constrainIndex(const int idx, const int N) const;

  const int maxNx;
  const int maxNy;
  const int maxNz;

  std::array<double, X * Y * Z> m_data;
};

#include "../src/GridData.cpp"

using GridDataScalar = GridData<Nx, Ny, Nz>;
using GridDataX = GridData<Nx + 1, Ny, Nz>;
using GridDataY = GridData<Nx, Ny + 1, Nz>;
using GridDataZ = GridData<Nx, Ny, Nz + 1>;