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
  int constrainIndex(int idx, int N);

private:
  double scalar[Nx * Ny * Nz];
  int maxNx;
  int maxNy;
  int maxNz;
};

class GridDataX : public GridData
{
public:
  GridDataX();
  virtual ~GridDataX();

  double &operator()(int i, int j, int k) override;

private:
  double mU[(Nx + 1) * Ny * Nz];
  int maxNx;
  int maxNy;
  int maxNz;
};

class GridDataY : public GridData
{
public:
  GridDataY();
  virtual ~GridDataY();

  double &operator()(int i, int j, int k) override;

private:
  double mV[Nx * (Ny + 1) * Nz];
  int maxNx;
  int maxNy;
  int maxNz;
};

class GridDataZ : public GridData
{
public:
  GridDataZ();
  virtual ~GridDataZ();

  double &operator()(int i, int j, int k) override;

private:
  double mW[Nx * Ny * (Nz + 1)];
  int maxNx;
  int maxNy;
  int maxNz;
};