#pragma once
#include "Vec3.hpp"
#include "constants.hpp"
#include "GridData.hpp"

#define FOR_EACH_CELL            \
  for (int k = 0; k < Nz; ++k)   \
    for (int j = 0; j < Ny; ++j) \
      for (int i = 0; i < Nx; ++i)

#define FOR_EACH_FACE_X          \
  for (int k = 0; k < Nz; ++k)   \
    for (int j = 0; j < Ny; ++j) \
      for (int i = 0; i < Nx + 1; ++i)

#define FOR_EACH_FACE_Y              \
  for (int k = 0; k < Nz; ++k)       \
    for (int j = 0; j < Ny + 1; ++j) \
      for (int i = 0; i < Nx; ++i)

#define FOR_EACH_FACE_Z            \
  for (int k = 0; k < Nz + 1; ++k) \
    for (int j = 0; j < Ny; ++j)   \
      for (int i = 0; i < Nx; ++i)

class MACGrid
{
public:
  MACGrid();
  ~MACGrid();

  Vec3 getCenter(int i, int j, int k);
  Vec3 getVelocity(const Vec3 &pos);

  double getVelocityX(const Vec3 &pos);
  double getVelocityY(const Vec3 &pos);
  double getVelocityZ(const Vec3 &pos);
  double getDensity(const Vec3 &pos);
  double getTemperature(const Vec3 &pos);
  double getPressure(const Vec3 &pos);

  GridDataX u, u0;
  GridDataY v, v0;
  GridDataZ w, w0;
  GridData density, temperature, pressure;
  double avg_u[SIZE], avg_v[SIZE], avg_w[SIZE];
  double omg_x[SIZE], omg_y[SIZE], omg_z[SIZE], vort[SIZE];
  double fx[SIZE], fy[SIZE], fz[SIZE];
};
