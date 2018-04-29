#include <algorithm>
#include <cassert>
#include "GridData.hpp"

/* GridData */
GridData::GridData() : scalar() {}
GridData::~GridData() {}
double &GridData::operator()(int i, int j, int k)
{
    assert((i >= 0 || i < Nx) || (j >= 0 || j < Ny) || (k >= 0 || k < Nz));
    return scalar[i + Nx * j + Nx * Ny * k];
}

double *GridData::getScalarPtr()
{
    return scalar;
}

double GridData::interp(const Vec3 &pt)
{
    Vec3 pos;
    pos[0] = std::min(std::max(0.0, pt[0] - VOXEL_SIZE * 0.5), (double)Nx);
    pos[1] = std::min(std::max(0.0, pt[1] - VOXEL_SIZE * 0.5), (double)Ny);
    pos[2] = std::min(std::max(0.0, pt[2] - VOXEL_SIZE * 0.5), (double)Nz);

    int i = (int)(pos[0] / VOXEL_SIZE);
    int j = (int)(pos[1] / VOXEL_SIZE);
    int k = (int)(pos[2] / VOXEL_SIZE);

    double scale = 1.0 / VOXEL_SIZE;
    double fractx = scale * (pos[0] - i * VOXEL_SIZE);
    double fracty = scale * (pos[1] - j * VOXEL_SIZE);
    double fractz = scale * (pos[2] - k * VOXEL_SIZE);

    assert(fractx < 1.0 && fractx >= 0);
    assert(fracty < 1.0 && fracty >= 0);
    assert(fractz < 1.0 && fractz >= 0);

    // Y @ low X, low Z:
    double tmp1 = (*this)(i, j, k);
    double tmp2 = (*this)(i, j + 1, k);
    // Y @ high X, low Z:
    double tmp3 = (*this)(i + 1, j, k);
    double tmp4 = (*this)(i + 1, j + 1, k);

    // Y @ low X, high Z:
    double tmp5 = (*this)(i, j, k + 1);
    double tmp6 = (*this)(i, j + 1, k + 1);
    // Y @ high X, high Z:
    double tmp7 = (*this)(i + 1, j, k + 1);
    double tmp8 = (*this)(i + 1, j + 1, k + 1);

    // Y @ low X, low Z
    double tmp12 = (1 - fracty) * tmp1 + fracty * tmp2;
    // Y @ high X, low Z
    double tmp34 = (1 - fracty) * tmp3 + fracty * tmp4;

    // Y @ low X, high Z
    double tmp56 = (1 - fracty) * tmp5 + fracty * tmp6;
    // Y @ high X, high Z
    double tmp78 = (1 - fracty) * tmp7 + fracty * tmp8;

    // X @ low Z
    double tmp1234 = (1 - fractx) * tmp12 + fractx * tmp34;
    // X @ high Z
    double tmp5678 = (1 - fractx) * tmp56 + fractx * tmp78;

    // Z
    double tmp = (1 - fractz) * tmp1234 + fractz * tmp5678;
    return tmp;
}

/* GridDataX */
GridDataX::GridDataX() : GridData(), mU() {}
GridDataX::~GridDataX() {}
double &GridDataX::operator()(int i, int j, int k)
{
    assert((i >= 0 || i <= Nx) || (j >= 0 || j < Ny) || (k >= 0 || k < Nz));
    return mU[i + (Nx + 1) * j + (Nx + 1) * Ny * k];
}

/* GridDataY */
GridDataY::GridDataY() : GridData(), mV() {}
GridDataY::~GridDataY() {}
double &GridDataY::operator()(int i, int j, int k)
{
    assert((i >= 0 || i < Nx) || (j >= 0 || j <= Ny) || (k >= 0 || k < Nz));
    return mV[i + Nx * j + Nx * (Ny + 1) * k];
}

/* GridDataZ */
GridDataZ::GridDataZ() : GridData(), mW() {}
GridDataZ::~GridDataZ() {}
double &GridDataZ::operator()(int i, int j, int k)
{
    assert((i >= 0 || i < Nx) || (j >= 0 || j < Ny) || (k >= 0 || k <= Nz));
    return mW[i + Nx * j + Nx * Ny * k];
}
