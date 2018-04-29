#pragma once
#include "constants.hpp"
#include "Vec3.hpp"

class GridData
{
public:
  GridData();
  ~GridData();

  virtual double &operator()(int i, int j, int k);
  virtual const double operator()(int i, int j, int k) const;
  virtual double *getScalarPtr();

  double interp(const Vec3 &pos);

private:
  double scalar[Nx * Ny * Nz];
};

class GridDataX : public GridData
{
public:
  GridDataX();
  ~GridDataX();

  double &operator()(int i, int j, int k) override;
  const double operator()(int i, int j, int k) const override;

private:
  double mU[(Nx + 1) * Ny * Nz];
};

class GridDataY : public GridData
{
public:
  GridDataY();
  ~GridDataY();

  double &operator()(int i, int j, int k) override;
  const double operator()(int i, int j, int k) const override;

private:
  double mV[Nx * (Ny + 1) * Nz];
};

class GridDataZ : public GridData
{
public:
  GridDataZ();
  ~GridDataZ();

  double &operator()(int i, int j, int k) override;
  const double operator()(int i, int j, int k) const override;

private:
  double mW[Nx * Ny * (Nz + 1)];
};