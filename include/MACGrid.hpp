#pragma once
#include "Vec3.hpp"
#include "constants.hpp"

#define FOR_EACH_CELL            \
  for (int k = 0; k < Nz; ++k)   \
    for (int j = 0; j < Ny; ++j) \
      for (int i = 0; i < Nx; ++i)

#define FOR_EACH_FACE                \
  for (int k = 0; k < Nz + 1; ++k)   \
    for (int j = 0; j < Ny + 1; ++j) \
      for (int i = 0; i < Nx + 1; ++i)

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

  Vec3 getVelocity(const Vec3 &pt);
  double getVelocityX(const Vec3 &pt);
  double getVelocityY(const Vec3 &pt);
  double getVelocityZ(const Vec3 &pt);
  Vec3 getCenter(int i, int j, int k);

  double u[MACSIZE_X], v[MACSIZE_Y], w[MACSIZE_Z];
  double u0[MACSIZE_X], v0[MACSIZE_Y], w0[MACSIZE_Z];

  double avg_u[SIZE], avg_v[SIZE], avg_w[SIZE];
  double omg_x[SIZE], omg_y[SIZE], omg_z[SIZE], omg_length[SIZE];
  double eta_x[SIZE], eta_y[SIZE], eta_z[SIZE];

  double density[SIZE], temperature[SIZE], pressure[SIZE];
  double fx[SIZE], fy[SIZE], fz[SIZE];
};
