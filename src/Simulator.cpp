#include <cmath>
#include <random>
#include "Simulator.hpp"

Simulator::Simulator(Voxels *voxels) : m_voxels(voxels), A(SIZE, SIZE), b(SIZE), x(SIZE)
{
    // nnz size is estimated by 7*SIZE because there are 7 nnz elements in a row.(center and neighbor 6)
    tripletList.reserve(7 * SIZE);

    /* add set temperature */
    std::random_device rnd;
    std::mt19937 mt(rnd());
    std::uniform_real_distribution<double> rand1(0, 1);
    for (int k = 0; k < N; ++k)
    {
        for (int i = 0; i < N; ++i)
        {
            m_voxels->temp[POS(i, N, k)] = 200 * rand1(mt);
        }
    }

    addSource();
    /* set emitter velocity */
    for (int k = N / 2 - SOURCE_SIZE / 2; k < N / 2 + SOURCE_SIZE / 2; ++k)
    {
        for (int j = 0; j < SOURCE_SIZE; ++j)
        {
            for (int i = N / 2 - SOURCE_SIZE / 2; i < N / 2 + SOURCE_SIZE / 2; ++i)
            {
                m_voxels->v[POSV(i, j, k)] = INIT_VELOCITY;
                m_voxels->v0[POSV(i, j, k)] = m_voxels->v[POSV(i, j, k)];
            }
        }
    }
}

Simulator::~Simulator()
{
}

void Simulator::update()
{
    //resetForce();
    averageVelocity();
    //calVorticity();
    addForce();
    advectVelocity();
    //calPressure();
    //applyPressureTerm();
    averageVelocity();
    advectScalar();
}

/* private */
void Simulator::addSource()
{
    for (int k = N / 2 - SOURCE_SIZE / 2; k < N / 2 + SOURCE_SIZE / 2; ++k)
    {
        for (int j = 0; j < SOURCE_SIZE; ++j)
        {
            for (int i = N / 2 - SOURCE_SIZE / 2; i < N / 2 + SOURCE_SIZE / 2; ++i)
            {
                m_voxels->dens[POS(i, j, k)] = INIT_DENSITY;
            }
        }
    }
}

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

void Simulator::averageVelocity()
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
}

void Simulator::calVorticity()
{
    for (int k = 0; k < N; ++k)
    {
        for (int j = 0; j < N; ++j)
        {
            for (int i = 0; i < N; ++i)
            {
                if (j > 0 && j < N - 1 && k > 0 || k < N - 1)
                {
                    m_voxels->omg_x[POS(i, j, k)] = (m_voxels->avg_w[POS(i, j + 1, k)] - m_voxels->avg_w[POS(i, j - 1, k)] - m_voxels->avg_v[POS(i, j, k + 1)] + m_voxels->avg_v[POS(i, j, k - 1)]) * 0.5 * N / (double)LENGTH;
                }
                if (k > 0 && k < N - 1 && i > 0 || i < N - 1)
                {
                    m_voxels->omg_y[POS(i, j, k)] = (m_voxels->avg_u[POS(i, j, k + 1)] - m_voxels->avg_u[POS(i, j, k - 1)] - m_voxels->avg_w[POS(i + 1, j, k)] + m_voxels->avg_w[POS(i - 1, j, k)]) * 0.5 * N / (double)LENGTH;
                }
                if (i > 0 && i < N - 1 && j > 0 || j < N - 1)
                {
                    m_voxels->omg_z[POS(i, j, k)] = (m_voxels->avg_v[POS(i + 1, j, k)] - m_voxels->avg_v[POS(i - 1, j, k)] - m_voxels->avg_u[POS(i, j + 1, k)] + m_voxels->avg_u[POS(i, j - 1, k)]) * 0.5 * N / (double)LENGTH;
                }
                m_voxels->omg_length[POS(i, j, k)] = l2norm(m_voxels->omg_x[POS(i, j, k)], m_voxels->omg_y[POS(i, j, k)], m_voxels->omg_z[POS(i, j, k)]);
            }
        }
    }

    for (int k = 0; k < N; ++k)
    {
        for (int j = 0; j < N; ++j)
        {
            for (int i = 0; i < N; ++i)
            {
                if (i < N - 1)
                {
                    m_voxels->eta_x[POS(i, j, k)] = (m_voxels->omg_length[POS(i + 1, j, k)] - m_voxels->omg_length[POS(i, j, k)]) * N / (double)LENGTH;
                }
                if (j < N - 1)
                {
                    m_voxels->eta_y[POS(i, j, k)] = (m_voxels->omg_length[POS(i, j + 1, k)] - m_voxels->omg_length[POS(i, j, k)]) * N / (double)LENGTH;
                }
                if (k < N - 1)
                {
                    m_voxels->eta_z[POS(i, j, k)] = (m_voxels->omg_length[POS(i, j, k + 1)] - m_voxels->omg_length[POS(i, j, k)]) * N / (double)LENGTH;
                }
                double norm = l2norm(m_voxels->eta_x[POS(i, j, k)], m_voxels->eta_y[POS(i, j, k)], m_voxels->eta_z[POS(i, j, k)]);
                if (norm != 0)
                {
                    m_voxels->eta_x[POS(i, j, k)] /= norm;
                    m_voxels->eta_y[POS(i, j, k)] /= norm;
                    m_voxels->eta_z[POS(i, j, k)] /= norm;
                }
            }
        }
    }

    for (int k = 0; k < N; ++k)
    {
        for (int j = 0; j < N; ++j)
        {
            for (int i = 0; i < N; ++i)
            {
                m_voxels->fx[POS(i, j, k)] += VORT_EPS * (N / (double)LENGTH) * (m_voxels->eta_y[POS(i, j, k)] * m_voxels->omg_z[POS(i, j, k)] - m_voxels->eta_z[POS(i, j, k)] * m_voxels->omg_y[POS(i, j, k)]);
                m_voxels->fy[POS(i, j, k)] += VORT_EPS * (N / (double)LENGTH) * (m_voxels->eta_z[POS(i, j, k)] * m_voxels->omg_x[POS(i, j, k)] - m_voxels->eta_x[POS(i, j, k)] * m_voxels->omg_z[POS(i, j, k)]);
                m_voxels->fz[POS(i, j, k)] += VORT_EPS * (N / (double)LENGTH) * (m_voxels->eta_x[POS(i, j, k)] * m_voxels->omg_y[POS(i, j, k)] - m_voxels->eta_y[POS(i, j, k)] * m_voxels->omg_x[POS(i, j, k)]);
            }
        }
    }
}

void Simulator::addForce()
{
    for (int k = 0; k < N; ++k)
    {
        for (int j = 0; j < N; ++j)
        {
            for (int i = 0; i < N; ++i)
            {
                if (i < N - 1)
                {
                    m_voxels->u[POSU(i + 1, j, k)] += DT * (m_voxels->fx[POS(i, j, k)] + m_voxels->fx[POS(i + 1, j, k)]) * 0.5;
                }
                if (j < N - 1)
                {
                    m_voxels->v[POSV(i, j + 1, k)] += DT * (m_voxels->fy[POS(i, j, k)] + m_voxels->fy[POS(i, j + 1, k)]) * 0.5;
                }
                if (k < N - 1)
                {
                    m_voxels->w[POSW(i, j, k + 1)] += DT * (m_voxels->fz[POS(i, j, k)] + m_voxels->fz[POS(i, j, k + 1)]) * 0.5;
                }

                m_voxels->u0[POSU(i + 1, j, k)] = m_voxels->u[POSU(i + 1, j, k)];
                m_voxels->v0[POSV(i, j + 1, k)] = m_voxels->v[POSV(i, j + 1, k)];
                m_voxels->w0[POSW(i, j, k + 1)] = m_voxels->w[POSW(i, j, k + 1)];
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
    tripletList.clear();
    A.setZero();
    b.setZero();
    x.setZero();

    double coeff = LENGTH * RHO / ((double)N * DT);

    for (int k = 0; k < N; ++k)
    {
        for (int j = 0; j < N; ++j)
        {
            for (int i = 0; i < N; ++i)
            {
                double F[6] = {k > 0, j > 0, i > 0, i < N, j < N, k < N};
                double D[6] = {-1.0, -1.0, -1.0, 1.0, 1.0, 1.0};
                double U[6];
                U[0] = (k > 0) ? m_voxels->w[POSW(i, j, k - 1)] : 0.0;
                U[1] = (j > 0) ? m_voxels->v[POSV(i, j - 1, k)] : 0.0;
                U[2] = (i > 0) ? m_voxels->u[POSU(i - 1, j, k)] : 0.0;
                U[3] = m_voxels->u[POSU(i + 1, j, k)];
                U[4] = m_voxels->v[POSV(i, j + 1, k)];
                U[5] = m_voxels->w[POSW(i, j, k + 1)];
                double sum_F = 0.0;

                for (int n = 0; n < 6; ++n)
                {
                    sum_F += F[n];
                    b(POS(i, j, k)) += D[n] * F[n] * U[n];
                }
                b(POS(i, j, k)) *= coeff;

                if (k > 0)
                {
                    tripletList.push_back(T(POS(i, j, k), POS(i, j, k - 1), F[0]));
                }
                if (j > 0)
                {
                    tripletList.push_back(T(POS(i, j, k), POS(i, j - 1, k), F[1]));
                }
                if (i > 0)
                {
                    tripletList.push_back(T(POS(i, j, k), POS(i - 1, j, k), F[2]));
                }

                tripletList.push_back(T(POS(i, j, k), POS(i, j, k), -sum_F));

                if (i < N - 1)
                {
                    tripletList.push_back(T(POS(i, j, k), POS(i + 1, j, k), F[3]));
                }
                if (j < N - 1)
                {
                    tripletList.push_back(T(POS(i, j, k), POS(i, j + 1, k), F[4]));
                }
                if (k < N - 1)
                {
                    tripletList.push_back(T(POS(i, j, k), POS(i, j, k + 1), F[5]));
                }
            }
        }
    }

    A.setFromTriplets(tripletList.begin(), tripletList.end());

    /* solve sparse lenear system by ICCG */
    ICCG.compute(A);
    if (ICCG.info() == Eigen::Success)
    {
        std::cout << "SUCCESS: Convergence" << std::endl;
    }
    else
    {
        std::cout << "FAILED: No Convergence" << std::endl;
    }
    x = ICCG.solve(b);
    std::cout << "#iterations:     " << ICCG.iterations() << std::endl;
    std::cout << "estimated error: " << ICCG.error() << std::endl;

    Eigen::Map<Eigen::VectorXd>(m_voxels->pressure, SIZE) = x;
}

void Simulator::applyPressureTerm()
{
    for (int k = 0; k < N; ++k)
    {
        for (int j = 0; j < N; ++j)
        {
            for (int i = 0; i < N; ++i)
            {
                if (i < N - 1)
                {
                    m_voxels->u[POSU(i, j, k)] -= DT * (m_voxels->pressure[POS(i + 1, j, k)] - m_voxels->pressure[POS(i, j, k)]) * N / (double)LENGTH;
                }
                if (j < N - 1)
                {
                    m_voxels->v[POSV(i, j, k)] -= DT * (m_voxels->pressure[POS(i, j + 1, k)] - m_voxels->pressure[POS(i, j, k)]) * N / (double)LENGTH;
                }
                if (k < N - 1)
                {
                    m_voxels->w[POSW(i, j, k)] -= DT * (m_voxels->pressure[POS(i, j, k + 1)] - m_voxels->pressure[POS(i, j, k)]) * N / (double)LENGTH;
                }
            }
        }
    }
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

                x = x - DT * interp(x, y, z, m_voxels->avg_u, N, N, N);
                y = y - DT * interp(x, y, z, m_voxels->avg_v, N, N, N);
                z = z - DT * interp(x, y, z, m_voxels->avg_w, N, N, N);

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

    double c[8] = {(1.0 - x) * (1.0 - y) * (1.0 - z), (1.0 - x) * (1.0 - y) * z, (1.0 - x) * y * (1.0 - z), x * (1.0 - y) * (1.0 - z), (1.0 - x) * y * z, x * (1.0 - y) * z, x * y * (1.0 - z), x * y * z};

    double ret = 0.0;
    for (int i = 0; i < 8; ++i)
    {
        ret += c[i] * f[i];
    }

    return ret;
}
