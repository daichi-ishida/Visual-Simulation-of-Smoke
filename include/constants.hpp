#pragma once
#include <cassert>

constexpr double VOXEL_SIZE = 1.0;
constexpr int Nx = 24, Ny = 48, Nz = 24;

constexpr int SOURCE_SIZE_X = 6;
constexpr int SOURCE_SIZE_Y = 2;
constexpr int SOURCE_SIZE_Z = 6;
constexpr int SOURCE_Y_MERGIN = 4;

constexpr double DT = 0.02;
constexpr double RHO = 1.0;
constexpr double INIT_DENSITY = 1.0;
constexpr double INIT_VELOCITY = 80.0;
constexpr double VORT_EPS = 0.55;
constexpr double ALPHA = 9.8;
constexpr double BETA = 12.0;
constexpr double T_AMP = 5.0;
constexpr double T_AMBIENT = 50.0;
constexpr double EMIT_DURATION = 2.0;
constexpr double FINISH_TIME = 6.0;

constexpr int SIZE = Nx * Ny * Nz;

constexpr int POS(int i, int j, int k)
{
    assert((i >= 0 || i < Nx) || (j >= 0 || j < Ny) || (k >= 0 || k < Nz));
    return i + Nx * j + Nx * Ny * k;
}