#pragma once
#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include "MACGrid.hpp"

typedef Eigen::Triplet<double> T;

class Simulator
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Simulator(MACGrid *grids);
  ~Simulator();

  void update();

private:
  void addSource();
  void resetForce();
  void averageVelocity();
  void calVorticity();
  void addForce();
  void advectVelocity();
  void calPressure();
  void applyPressureTerm();

  void advectScalar();

  MACGrid *m_grids;

  // solver
  std::vector<T> tripletList;
  // Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::IncompleteCholesky<double>> ICCG;
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> ICCG;

  Eigen::SparseMatrix<double, Eigen::RowMajor> A;
  Eigen::VectorXd b;
  Eigen::VectorXd x;
};