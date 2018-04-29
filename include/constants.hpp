#pragma once
#include <cassert>

constexpr double VOXEL_SIZE = 1.0;
constexpr int Nx = 32, Ny = 32, Nz = 32;

constexpr int SOURCE_SIZE_X = Nx / 3;
constexpr int SOURCE_SIZE_Y = Nx / 3;
constexpr int SOURCE_SIZE_Z = Nx / 3;

constexpr double DT = 0.1;
constexpr double RHO = 1.0;
constexpr double INIT_DENSITY = 1.0;
constexpr double INIT_VELOCITY = 2.0;
constexpr double VORT_EPS = 0.6;
constexpr double GRAVITY_Y = 9.8;
constexpr double T_AMBIENT = 30.0;
constexpr double FINISH_TIME = 20.0;

constexpr int SIZE = Nx * Ny * Nz;

constexpr int POS(int i, int j, int k)
{
    assert((i >= 0 || i < Nx) || (j >= 0 || j < Ny) || (k >= 0 || k < Nz));
    return i + Nx * j + Nx * Ny * k;
}