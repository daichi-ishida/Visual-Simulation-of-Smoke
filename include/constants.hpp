#pragma once
#include <iostream>
#include <string>

/* Scene Constants */
static const char *TITLE = "Visual Simulation of Smoke";
constexpr int LENGTH = 640;

/* Simulator Constants */
constexpr int N = 32;
constexpr double DT = 0.01;
constexpr double VISCOSITY = 0.001;
constexpr double VORT_EPS = 10.0;
constexpr double GRAVITY_Y = 9.8;
constexpr double T_AMBIENT = 30.0;
constexpr double FINISH_TIME = 3.0;

constexpr int SIZE = N * N * N;
constexpr int MAC_SIZE = N * N * (N + 1);

constexpr int POS(int i, int j, int k) { return i + N * j + N * N * k; };

constexpr int POSU(int i, int j, int k) { return i + (N + 1) * j + N * (N + 1) * k; };
constexpr int POSV(int i, int j, int k) { return i + N * j + (N + 1) * N * k; };
constexpr int POSW(int i, int j, int k) { return i + N * j + N * N * k; };

constexpr double l2norm(double x, double y, double z) { return std::sqrt(x * x + y * y + z * z); };

enum ETag
{
    E_FLUID = 0,
    E_BOUNDARY = 1,
    E_INSIDE_OBJ = 2
}