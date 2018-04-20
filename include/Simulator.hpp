#pragma once
#include <fftw3.h>

#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>
#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>

#include "Voxels.hpp"

class Simulator
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Simulator(Voxels *voxels);
  ~Simulator();

  void update();

private:
  void addSource();
  void resetForce();
  void calVorticity();
  void addForce();
  void advectVelocity();
  void calPressure();
  void applyPressureTerm();

  void advectScalar();

  double interp(double x, double y, double z, double q[], unsigned int Nx, unsigned int Ny, unsigned int Nz);

  Voxels *m_voxels;

  // solver
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::IncompleteCholesky<double>> ICCG;
  Eigen::SparseMatrix<double, Eigen::RowMajor> A;
  Eigen::VectorXd b;
  Eigen::VectorXd x;
};