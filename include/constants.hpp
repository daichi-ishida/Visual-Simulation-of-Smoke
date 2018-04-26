#pragma once
#include <cmath>
#include <iostream>
#include <string>

/* Scene Constants */
constexpr int LENGTH = 640;

/* Simulator Constants */
constexpr int N = 32;
constexpr int SOURCE_SIZE = N / 3;
constexpr int SOURCE_MARGIN = N / 10;
constexpr double DT = 0.1;
constexpr double RHO = 1.0;
constexpr double INIT_DENSITY = 1.0;
constexpr double INIT_VELOCITY = 2.0;
constexpr double VORT_EPS = 0.01;
constexpr double GRAVITY_Y = 9.8;
constexpr double T_AMBIENT = 30.0;
constexpr double FINISH_TIME = 20.0;

constexpr int SIZE = N * N * N;
constexpr int MAC_SIZE = N * N * (N + 1);

constexpr int POS(int i, int j, int k) { return i + N * j + N * N * k; };

constexpr int POSU(int i, int j, int k) { return i + (N + 1) * j + N * (N + 1) * k; };
constexpr int POSV(int i, int j, int k) { return i + N * j + (N + 1) * N * k; };
constexpr int POSW(int i, int j, int k) { return i + N * j + N * N * k; };

constexpr double l2norm(double x, double y, double z) { return std::sqrt(x * x + y * y + z * z); };

enum EMode
{
    E_U = 0,
    E_V = 1,
    E_W = 2
};