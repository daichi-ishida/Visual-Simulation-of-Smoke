#pragma once
#include <cassert>

constexpr double VOXEL_SIZE = 1.0;
constexpr int Nx = 32, Ny = 64, Nz = 32;

constexpr int SOURCE_SIZE_X = Nx / 3;
constexpr int SOURCE_SIZE_Y = Nx / 3;
constexpr int SOURCE_SIZE_Z = Nx / 3;

constexpr double DT = 0.1;
constexpr double RHO = 1.0;
constexpr double INIT_DENSITY = 1.0;
constexpr double INIT_VELOCITY = 2.0;
constexpr double VORT_EPS = 0.01;
constexpr double GRAVITY_Y = 9.8;
constexpr double T_AMBIENT = 30.0;
constexpr double FINISH_TIME = 20.0;

constexpr int SIZE = Nx * Ny * Nz;
constexpr int MACSIZE_X = (Nx + 1) * Ny * Nz;
constexpr int MACSIZE_Y = Nx * (Ny + 1) * Nz;
constexpr int MACSIZE_Z = Nx * Ny * (Nz + 1);

constexpr int POS(int i, int j, int k)
{
    assert((i >= 0 || i < Nx) || (j >= 0 || j < Ny) || (k >= 0 || k < Nz));
    return i + Nx * j + Nx * Ny * k;
}

constexpr int POSU(int i, int j, int k)
{
    assert((i >= 0 || i <= Nx) || (j >= 0 || j < Ny) || (k >= 0 || k < Nz));
    return i + (Nx + 1) * j + (Nx + 1) * Ny * k;
}
constexpr int POSV(int i, int j, int k)
{
    assert((i >= 0 || i < Nx) || (j >= 0 || j <= Ny) || (k >= 0 || k < Nz));
    return i + Nx * j + Nx * (Ny + 1) * k;
}
constexpr int POSW(int i, int j, int k)
{
    assert((i >= 0 || i < Nx) || (j >= 0 || j < Ny) || (k >= 0 || k <= Nz));
    return i + Nx * j + Nx * Ny * k;
}