#pragma once
#include "constants.hpp"
#include "Vec3.hpp"

class GridData
{
public:
  GridData();
  ~GridData();

  virtual double &operator()(int i, int j, int k);
  virtual double *getScalarPtr();

  double interp(const Vec3 &pt);

protected:
  double linearInterpolation(const Vec3 &pt);
  double monotonicCubicInterpolation(const Vec3 &pt);
  double axis_monotonicCubicInterpolation(double f[], double fract);
  int sign(double a);

private:
  double scalar[Nx * Ny * Nz];
};

class GridDataX : public GridData
{
public:
  GridDataX();
  virtual ~GridDataX();

  virtual double &operator()(int i, int j, int k);

private:
  double mU[(Nx + 1) * Ny * Nz];
};

class GridDataY : public GridData
{
public:
  GridDataY();
  virtual ~GridDataY();

  virtual double &operator()(int i, int j, int k);

private:
  double mV[Nx * (Ny + 1) * Nz];
};

class GridDataZ : public GridData
{
public:
  GridDataZ();
  virtual ~GridDataZ();

  virtual double &operator()(int i, int j, int k);

private:
  double mW[Nx * Ny * (Nz + 1)];
};