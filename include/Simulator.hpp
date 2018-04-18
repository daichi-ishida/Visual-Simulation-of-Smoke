#pragma once
#include <fftw3.h>

#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>
#pragma once
#include "Voxels.hpp"

class Simulator
{
public:
  Simulator(Voxels *voxels);
  ~Simulator();

  void update();

private:
  void resetForce();
  void calVorticity();
  void addForce();
  void advectVelocity();
  void calPressure();
  void calPressureGradient();
  void applyPressureTerm();

  void advectScalar();

  double interp(double x, double y, double z, double q[], unsigned int Nx, unsigned int Ny, unsigned int Nz);

  Voxels *m_voxels;
};