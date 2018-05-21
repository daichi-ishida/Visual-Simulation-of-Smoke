#pragma once
#include <cassert>

enum E_METHOD
{
    E_LINEAR = 0,
    E_MONOTONIC_CUBIC = 1
};

/* Simulation Constants */
constexpr double VOXEL_SIZE = 1.0;
constexpr int Nx = 50, Ny = 50, Nz = 50;
constexpr E_METHOD INTERPOLATION_METHOD = E_MONOTONIC_CUBIC;

constexpr int SOURCE_SIZE_X = 6;
constexpr int SOURCE_SIZE_Y = 2;
constexpr int SOURCE_SIZE_Z = 6;
constexpr int SOURCE_Y_MERGIN = 4;

constexpr double DT = 0.02;
constexpr double INIT_DENSITY = 1.0;
constexpr double INIT_VELOCITY = 60.0;
constexpr double VORT_EPS = 0.25;
constexpr double ALPHA = 9.8;
constexpr double BETA = 15.0;
constexpr double T_AMP = 5.0;
constexpr double T_AMBIENT = 50.0;
constexpr double EMIT_DURATION = 1.0;
constexpr double FINISH_TIME = 6.0;

constexpr int SIZE = Nx * Ny * Nz;

constexpr int POS(int i, int j, int k)
{
    assert((i >= 0 || i < Nx) || (j >= 0 || j < Ny) || (k >= 0 || k < Nz));
    return i + Nx * j + Nx * Ny * k;
}

/* Scene Constants */
constexpr int WIN_WIDTH = 500;
constexpr int WIN_HEIGHT = 500;
static const char *WIN_TITLE = "Test";
constexpr float RADIUS_THRESHOLD = 0.5f;

#ifdef _OPENMP
#include <omp.h>
#define OPENMP_FOR _Pragma("omp parallel for")
#define OPENMP_SECTION _Pragma("omp section")
#define OPENMP_BEGIN        \
    _Pragma("omp parallel") \
    {
#define OPENMP_END }
#define OPENMP_FOR_P _Pragma("omp for")
#else
#define OPENMP_FOR
#define OPENMP_SECTION
#define OPENMP_BEGIN
#define OPENMP_END
#define OPENMP_FOR_P
#endif