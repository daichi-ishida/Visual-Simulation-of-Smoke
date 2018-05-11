#include <algorithm>
#include <cassert>
#include <cmath>
#include "GridData.hpp"

/* GridData */
GridData::GridData() : scalar(), maxNx(Nx - 1), maxNy(Ny - 1), maxNz(Nz - 1) {}
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
    switch (INTERPOLATION_METHOD)
    {
    case E_LINEAR:
        return linearInterpolation(pt);
    case E_MONOTONIC_CUBIC:
        return monotonicCubicInterpolation(pt);
    }
}

double GridData::linearInterpolation(const Vec3 &pt)
{
    Vec3 pos;
    // clamp position
    pos[0] = std::min(std::max(0.0, pt[0]), (double)maxNx * VOXEL_SIZE);
    pos[1] = std::min(std::max(0.0, pt[1]), (double)maxNy * VOXEL_SIZE);
    pos[2] = std::min(std::max(0.0, pt[2]), (double)maxNz * VOXEL_SIZE);

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

double GridData::monotonicCubicInterpolation(const Vec3 &pt)
{
    Vec3 pos;
    // clamp position
    pos[0] = std::min(std::max(0.0, pt[0]), (double)(maxNx - 1) * VOXEL_SIZE - 1e-6);
    pos[1] = std::min(std::max(0.0, pt[1]), (double)(maxNy - 1) * VOXEL_SIZE - 1e-6);
    pos[2] = std::min(std::max(0.0, pt[2]), (double)(maxNz - 1) * VOXEL_SIZE - 1e-6);

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

    // Z:
    double arr_z[4];
    for (int z = 0; z < 4; ++z)
    {
        // X @ Z_(z-1):
        double arr_x[4];
        for (int x = 0; x < 4; ++x)
        {
            // Y @ X_(x-1), Z_(z-1):
            int i1 = constrainIndex(i + x - 1, maxNx);
            int j1 = constrainIndex(j - 1, maxNy);
            int j2 = constrainIndex(j + 2, maxNy);
            int k1 = constrainIndex(k + z - 1, maxNz);

            double arr_y[4] = {(*this)(i1, j1, k1), (*this)(i1, j, k1), (*this)(i1, j + 1, k1), (*this)(i1, j2, k1)};
            arr_x[x] = axis_monotonicCubicInterpolation(arr_y, fracty);
        }
        arr_z[z] = axis_monotonicCubicInterpolation(arr_x, fractx);
    }

    return axis_monotonicCubicInterpolation(arr_z, fractz);
}

double GridData::axis_monotonicCubicInterpolation(double f[], double t)
{
    // f[0]:f_k-1, f[1]:f_k, f[2]:f_k+1, f[3]:f_k+2
    double delta = f[2] - f[1];
    double d0 = 0.5 * (f[2] - f[0]);
    double d1 = 0.5 * (f[3] - f[1]);

    // neccessary condition for monotonic
    d0 = sign(delta) * std::fabs(d0);
    d1 = sign(delta) * std::fabs(d1);

    double a0 = f[1];
    double a1 = d0;
    double a2 = 3 * delta - 2 * d0 - d1;
    double a3 = d0 + d1 - 2 * delta;
    return a3 * t * t * t + a2 * t * t + a1 * t + a0;
}

int GridData::sign(double a)
{
    // if a is positive, return 1
    // if a is negative, return -1
    // if a is zero, return 0
    return (a > 0) - (a < 0);
}

int GridData::constrainIndex(int idx, int N)
{
    if (idx == -1)
    {
        return 0;
    }
    if (idx > N)
    {
        return N;
    }
    return idx;
}

/* GridDataX */
GridDataX::GridDataX() : GridData(), maxNx(Nx), maxNy(Ny - 1), maxNz(Nz - 1), mU() {}
GridDataX::~GridDataX() {}
double &GridDataX::operator()(int i, int j, int k)
{
    assert((i >= 0 || i <= Nx) || (j >= 0 || j < Ny) || (k >= 0 || k < Nz));
    return mU[i + (Nx + 1) * j + (Nx + 1) * Ny * k];
}

/* GridDataY */
GridDataY::GridDataY() : GridData(), maxNx(Nx - 1), maxNy(Ny), maxNz(Nz - 1), mV() {}
GridDataY::~GridDataY() {}
double &GridDataY::operator()(int i, int j, int k)
{
    assert((i >= 0 || i < Nx) || (j >= 0 || j <= Ny) || (k >= 0 || k < Nz));
    return mV[i + Nx * j + Nx * (Ny + 1) * k];
}

/* GridDataZ */
GridDataZ::GridDataZ() : GridData(), maxNx(Nx - 1), maxNy(Ny - 1), maxNz(Nz), mW() {}
GridDataZ::~GridDataZ() {}
double &GridDataZ::operator()(int i, int j, int k)
{
    assert((i >= 0 || i < Nx) || (j >= 0 || j < Ny) || (k >= 0 || k <= Nz));
    return mW[i + Nx * j + Nx * Ny * k];
}
