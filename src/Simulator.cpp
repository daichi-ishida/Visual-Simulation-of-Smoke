#include <cmath>
#include <Eigen/Core>

#include "Simulator.hpp"

Simulator::Simulator(Voxels *voxels) : m_voxels(voxels)
{
}

Simulator::~Simulator()
{
}

void Simulator::update()
{
    resetForce();
    calVorticity();
    addForce();
    advectVelocity();
    calPressure();
    calPressureGradient();
    applyPressureTerm();

    advectScalar();
}

/* private */
void Simulator::resetForce()
{
    for (int k = 0; k < N; ++k)
    {
        for (int j = 0; j < N; ++j)
        {
            for (int i = 0; i < N; ++i)
            {
                m_voxels->fx[POS(i, j, k)] = 0.0;
                m_voxels->fy[POS(i, j, k)] = GRAVITY_Y * m_voxels->dens[POS(i, j, k)] - (m_voxels->temp[POS(i, j, k)] - T_AMBIENT);
                m_voxels->fz[POS(i, j, k)] = 0.0;
            }
        }
    }
}

void Simulator::calVorticity()
{
    for (int k = 0; k < N; ++k)
    {
        for (int j = 0; j < N; ++j)
        {
            for (int i = 0; i < N; ++i)
            {
                m_voxels->avg_u[POS(i, j, k)] = (m_voxels->u[POSU(i, j, k)] + m_voxels->u[POSU(i + 1, j, k)]) * 0.5;
                m_voxels->avg_v[POS(i, j, k)] = (m_voxels->v[POSV(i, j, k)] + m_voxels->v[POSV(i, j + 1, k)]) * 0.5;
                m_voxels->avg_w[POS(i, j, k)] = (m_voxels->w[POSW(i, j, k)] + m_voxels->w[POSW(i, j, k + 1)]) * 0.5;
            }
        }
    }

    for (int k = 1; k < N - 1; ++k)
    {
        for (int j = 1; j < N - 1; ++j)
        {
            for (int i = 1; i < N - 1; ++i)
            {
                m_voxels->omg_x[POS(i, j, k)] = (m_voxels->avg_w[POS(i, j + 1, k)] - m_voxels->avg_w[POS(i, j - 1, k)] - m_voxels->avg_v[POS(i, j, k + 1)] + m_voxels->avg_v[POS(i, j, k - 1)]) * 0.5 * N / LENGTH;
                m_voxels->omg_y[POS(i, j, k)] = (m_voxels->avg_u[POS(i, j, k + 1)] - m_voxels->avg_u[POS(i, j, k - 1)] - m_voxels->avg_w[POS(i + 1, j, k)] + m_voxels->avg_w[POS(i - 1, j, k)]) * 0.5 * N / LENGTH;
                m_voxels->omg_z[POS(i, j, k)] = (m_voxels->avg_v[POS(i + 1, j, k)] - m_voxels->avg_v[POS(i - 1, j, k)] - m_voxels->avg_u[POS(i, j + 1, k)] + m_voxels->avg_u[POS(i, j - 1, k)]) * 0.5 * N / LENGTH;
                m_voxels->omg_length[POS(i, j, k)] = l2norm(m_voxels->omg_x[POS(i, j, k)], m_voxels->omg_y[POS(i, j, k)], m_voxels->omg_z[POS(i, j, k)]);
            }
        }
    }

    for (int k = 0; k < N - 1; ++k)
    {
        for (int j = 0; j < N - 1; ++j)
        {
            for (int i = 0; i < N - 1; ++i)
            {
                m_voxels->eta_x[POS(i, j, k)] = (m_voxels->omg_length[POS(i + 1, j, k)] - m_voxels->omg_length[POS(i, j, k)]) * N / LENGTH;
                m_voxels->eta_y[POS(i, j, k)] = (m_voxels->omg_length[POS(i, j + 1, k)] - m_voxels->omg_length[POS(i, j, k)]) * N / LENGTH;
                m_voxels->eta_z[POS(i, j, k)] = (m_voxels->omg_length[POS(i, j, k + 1)] - m_voxels->omg_length[POS(i, j, k)]) * N / LENGTH;
                double norm = l2norm(m_voxels->eta_x[POS(i, j, k)], m_voxels->eta_y[POS(i, j, k)], m_voxels->eta_z[POS(i, j, k)]);
                m_voxels->eta_x[POS(i, j, k)] /= norm;
                m_voxels->eta_y[POS(i, j, k)] /= norm;
                m_voxels->eta_z[POS(i, j, k)] /= norm;
            }
        }
    }

    for (int k = 0; k < N - 1; ++k)
    {
        for (int j = 0; j < N - 1; ++j)
        {
            for (int i = 0; i < N - 1; ++i)
            {
                m_voxels->fx[POS(i, j, k)] += VORT_EPS * (N / LENGTH) * (m_voxels->eta_y[POS(i, j, k)] * m_voxels->omg_z[POS(i, j, k)] - m_voxels->eta_z[POS(i, j, k)] * m_voxels->omg_y[POS(i, j, k)]);
                m_voxels->fy[POS(i, j, k)] += VORT_EPS * (N / LENGTH) * (m_voxels->eta_z[POS(i, j, k)] * m_voxels->omg_x[POS(i, j, k)] - m_voxels->eta_x[POS(i, j, k)] * m_voxels->omg_z[POS(i, j, k)]);
                m_voxels->fz[POS(i, j, k)] += VORT_EPS * (N / LENGTH) * (m_voxels->eta_x[POS(i, j, k)] * m_voxels->omg_y[POS(i, j, k)] - m_voxels->eta_y[POS(i, j, k)] * m_voxels->omg_x[POS(i, j, k)]);
            }
        }
    }
}

void Simulator::addForce()
{
    for (int k = 0; k < N - 1; ++k)
    {
        for (int j = 0; j < N - 1; ++j)
        {
            for (int i = 0; i < N - 1; ++i)
            {
                m_voxels->u[POSU(i + 1, j, k)] += DT * (m_voxels->fx[POS(i, j, k)] + m_voxels->fx[POS(i + 1, j, k)]) * 0.5;
                m_voxels->v[POSV(i, j + 1, k)] += DT * (m_voxels->fy[POS(i, j, k)] + m_voxels->fy[POS(i, j + 1, k)]) * 0.5;
                m_voxels->w[POSW(i, j, k + 1)] += DT * (m_voxels->fz[POS(i, j, k)] + m_voxels->fz[POS(i, j, k + 1)]) * 0.5;
            }
        }
    }
}

void Simulator::advectVelocity()
{
    for (int k = 0; k < N; ++k)
    {
        for (int j = 0; j < N; ++j)
        {
            for (int i = 1; i < N; ++i)
            {
                double x = i * LENGTH / (double)N;
                double y = (j + 0.5) * LENGTH / (double)N;
                double z = (k + 0.5) * LENGTH / (double)N;

                x = x - DT * interp(x, y - 0.5 * LENGTH / (double)N, z - 0.5 * LENGTH / (double)N, m_voxels->u0, N + 1, N, N);
                y = y - DT * interp(x - 0.5 * LENGTH / (double)N, y, z - 0.5 * LENGTH / (double)N, m_voxels->v0, N, N + 1, N);
                z = z - DT * interp(x - 0.5 * LENGTH / (double)N, y - 0.5 * LENGTH / (double)N, z, m_voxels->w0, N, N, N + 1);

                m_voxels->u[POSU(i, j, k)] = interp(x, y - 0.5 * LENGTH / (double)N, z - 0.5 * LENGTH / (double)N, m_voxels->u0, N + 1, N, N);
            }
        }
    }
    for (int k = 0; k < N; ++k)
    {
        for (int j = 1; j < N; ++j)
        {
            for (int i = 0; i < N; ++i)
            {
                double x = (i + 0.5) * LENGTH / (double)N;
                double y = j * LENGTH / (double)N;
                double z = (k + 0.5) * LENGTH / (double)N;

                x = x - DT * interp(x, y - 0.5 * LENGTH / (double)N, z - 0.5 * LENGTH / (double)N, m_voxels->u0, N + 1, N, N);
                y = y - DT * interp(x - 0.5 * LENGTH / (double)N, y, z - 0.5 * LENGTH / (double)N, m_voxels->v0, N, N + 1, N);
                z = z - DT * interp(x - 0.5 * LENGTH / (double)N, y - 0.5 * LENGTH / (double)N, z, m_voxels->w0, N, N, N + 1);

                m_voxels->v[POSV(i, j, k)] = interp(x - 0.5 * LENGTH / (double)N, y, z - 0.5 * LENGTH / (double)N, m_voxels->v0, N, N + 1, N);
            }
        }
    }
    for (int k = 1; k < N; ++k)
    {
        for (int j = 0; j < N; ++j)
        {
            for (int i = 0; i < N; ++i)
            {
                double x = (i + 0.5) * LENGTH / (double)N;
                double y = (j + 0.5) * LENGTH / (double)N;
                double z = k * LENGTH / (double)N;

                x = x - DT * interp(x, y - 0.5 * LENGTH / (double)N, z - 0.5 * LENGTH / (double)N, m_voxels->u0, N + 1, N, N);
                y = y - DT * interp(x - 0.5 * LENGTH / (double)N, y, z - 0.5 * LENGTH / (double)N, m_voxels->v0, N, N + 1, N);
                z = z - DT * interp(x - 0.5 * LENGTH / (double)N, y - 0.5 * LENGTH / (double)N, z, m_voxels->w0, N, N, N + 1);

                m_voxels->w[POSW(i, j, k)] = interp(x - 0.5 * LENGTH / (double)N, y - 0.5 * LENGTH / (double)N, z, m_voxels->w0, N, N, N + 1);
            }
        }
    }
}

void Simulator::calPressure()
{
}

void Simulator::calPressureGradient()
{
}

void Simulator::applyPressureTerm()
{
}

void Simulator::advectScalar()
{
    for (unsigned int k = 0; k < N; ++k)
    {
        for (unsigned int j = 0; j < N; ++j)
        {
            for (unsigned int i = 0; i < N; ++i)
            {
                double x = i * LENGTH / (double)N;
                double y = j * LENGTH / (double)N;
                double z = k * LENGTH / (double)N;

                x = x - DT * interp(x, y, z, m_voxels->u, N, N, N);
                y = y - DT * interp(x, y, z, m_voxels->v, N, N, N);
                z = z - DT * interp(x, y, z, m_voxels->w, N, N, N);

                m_voxels->dens[POS(i, j, k)] = interp(x, y, z, m_voxels->dens, N, N, N);
                m_voxels->temp[POS(i, j, k)] = interp(x, y, z, m_voxels->temp, N, N, N);
            }
        }
    }
}

double Simulator::interp(double x, double y, double z, double q[], unsigned int Nx, unsigned int Ny, unsigned int Nz)
{
    x = std::fmax(0.0, std::fmin(Nx - 1 - 1e-6, N * x / (double)LENGTH));
    y = std::fmax(0.0, std::fmin(Ny - 1 - 1e-6, N * y / (double)LENGTH));
    z = std::fmax(0.0, std::fmin(Nz - 1 - 1e-6, N * z / (double)LENGTH));

    unsigned int i = x;
    unsigned int j = y;
    unsigned int k = z;

    double f[8] = {q[POS(i, j, k)], q[POS(i, j, k + 1)], q[POS(i, j + 1, k)], q[POS(i + 1, j, k)], q[POS(i, j + 1, k + 1)], q[POS(i + 1, j, k + 1)], q[POS(i + 1, j + 1, k)], q[POS(i + 1, j + 1, k + 1)]};

    x = x - i;
    y = y - j;
    z = z - k;

    double c[8] = {(1.0f - x) * (1.0f - y) * (1.0f - z), (1.0f - x) * (1.0f - y) * z, (1.0f - x) * y * (1.0f - z), x * (1.0f - y) * (1.0f - z), (1.0f - x) * y * z, x * (1.0f - y) * z, x * y * (1.0f - z), x * y * z};

    double ret = 0.0;
    for (int i = 0; i < 8; ++i)
    {
        ret += c[i] * f[i];
    }

    return ret;
}
