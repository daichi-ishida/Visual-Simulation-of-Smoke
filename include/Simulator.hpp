#pragma once
#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include "Voxels.hpp"

typedef Eigen::Triplet<double> T;

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
  std::vector<T> tripletList;
  // Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::IncompleteCholesky<double>> ICCG;
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> ICCG;

  Eigen::SparseMatrix<double, Eigen::RowMajor> A;
  Eigen::VectorXd b;
  Eigen::VectorXd x;
};