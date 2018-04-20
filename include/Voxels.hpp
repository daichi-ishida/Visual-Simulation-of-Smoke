#pragma once
#include <Eigen/Core>

#include "constants.hpp"

class Voxels
{
public:
  Voxels();
  ~Voxels();

  double u[MAC_SIZE], v[MAC_SIZE], w[MAC_SIZE];
  double u0[MAC_SIZE], v0[MAC_SIZE], w0[MAC_SIZE];

  double avg_u[SIZE], avg_v[SIZE], avg_w[SIZE];
  double omg_x[SIZE], omg_y[SIZE], omg_z[SIZE], omg_length[SIZE];
  double eta_x[SIZE], eta_y[SIZE], eta_z[SIZE];

  double dens[SIZE], temp[SIZE], pressure[SIZE];
  double fx[SIZE], fy[SIZE], fz[SIZE];
  bool is_fluid[SIZE];
};