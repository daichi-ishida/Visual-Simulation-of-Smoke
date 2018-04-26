#include <cmath>
#include <random>
#include "Simulator.hpp"

Simulator::Simulator(MACGrid *grids) : m_grids(grids), A(SIZE, SIZE), b(SIZE), x(SIZE)
{
    // nnz size is estimated by 7*SIZE because there are 7 nnz elements in a row.(center and neighbor 6)
    tripletList.reserve(7 * SIZE);

    /* add set temperatureerature */
    std::random_device rnd;
    std::mt19937 mt(rnd());
    std::uniform_real_distribution<double> rand1(0, 1);
    for (int k = 0; k < Nz; ++k)
    {
        for (int i = 0; i < Nx; ++i)
        {
            m_grids->temperature[POS(i, Ny, k)] = 200 * rand1(mt);
        }
    }

    addSource();
    /* set emitter velocity */
    for (int k = Nz / 2 - SOURCE_SIZE_Z / 2; k < Nz / 2 + SOURCE_SIZE_Z / 2; ++k)
    {
        for (int j = 0; j < SOURCE_SIZE_Y; ++j)
        {
            for (int i = Nx / 2 - SOURCE_SIZE_X / 2; i < Nx / 2 + SOURCE_SIZE_X / 2; ++i)
            {
                m_grids->v[POSV(i, j, k)] = INIT_VELOCITY;
                m_grids->v0[POSV(i, j, k)] = m_grids->v[POSV(i, j, k)];
            }
        }
    }
}

Simulator::~Simulator()
{
}

void Simulator::update()
{
    resetForce();
    averageVelocity();
    calVorticity();
    addForce();
    advectVelocity();
    calPressure();
    applyPressureTerm();
    averageVelocity();
    advectScalar();
}

/* private */
void Simulator::addSource()
{
    for (int k = Nz / 2 - SOURCE_SIZE_Z / 2; k < Nz / 2 + SOURCE_SIZE_Z / 2; ++k)
    {
        for (int j = 0; j < SOURCE_SIZE_Y; ++j)
        {
            for (int i = Nx / 2 - SOURCE_SIZE_X / 2; i < Nx / 2 + SOURCE_SIZE_X / 2; ++i)
            {
                m_grids->density[POS(i, j, k)] = INIT_DENSITY;
            }
        }
    }
}

void Simulator::resetForce()
{
    FOR_EACH_CELL
    {
        m_grids->fx[POS(i, j, k)] = 0.0;
        m_grids->fy[POS(i, j, k)] = GRAVITY_Y * m_grids->density[POS(i, j, k)] - (m_grids->temperature[POS(i, j, k)] - T_AMBIENT);
        m_grids->fz[POS(i, j, k)] = 0.0;
    }
}

void Simulator::averageVelocity()
{
    FOR_EACH_CELL
    {
        m_grids->avg_u[POS(i, j, k)] = (m_grids->u[POSU(i, j, k)] + m_grids->u[POSU(i + 1, j, k)]) * 0.5;
        m_grids->avg_v[POS(i, j, k)] = (m_grids->v[POSV(i, j, k)] + m_grids->v[POSV(i, j + 1, k)]) * 0.5;
        m_grids->avg_w[POS(i, j, k)] = (m_grids->w[POSW(i, j, k)] + m_grids->w[POSW(i, j, k + 1)]) * 0.5;
    }
}

void Simulator::calVorticity()
{
    FOR_EACH_CELL
    {
        // ignore boundary cells
        if (i == 0 || j == 0 || k == 0)
        {
            continue;
        }
        if (i == Nx - 1 || j == Ny - 1 || k == Nz - 1)
        {
            continue;
        }

        m_grids->omg_x[POS(i, j, k)] = (m_grids->avg_w[POS(i, j + 1, k)] - m_grids->avg_w[POS(i, j - 1, k)] - m_grids->avg_v[POS(i, j, k + 1)] + m_grids->avg_v[POS(i, j, k - 1)]) * 0.5 / VOXEL_SIZE;
        m_grids->omg_y[POS(i, j, k)] = (m_grids->avg_u[POS(i, j, k + 1)] - m_grids->avg_u[POS(i, j, k - 1)] - m_grids->avg_w[POS(i + 1, j, k)] + m_grids->avg_w[POS(i - 1, j, k)]) * 0.5 / VOXEL_SIZE;
        m_grids->omg_z[POS(i, j, k)] = (m_grids->avg_v[POS(i + 1, j, k)] - m_grids->avg_v[POS(i - 1, j, k)] - m_grids->avg_u[POS(i, j + 1, k)] + m_grids->avg_u[POS(i, j - 1, k)]) * 0.5 / VOXEL_SIZE;
    }

    FOR_EACH_CELL
    {
        // ignore boundary cells
        if (i == 0 || j == 0 || k == 0)
        {
            continue;
        }
        if (i == Nx - 1 || j == Ny - 1 || k == Nz - 1)
        {
            continue;
        }
        // compute gradient of vorticity using central differences
        double p, q;
        p = Vec3(m_grids->omg_x[POS(i + 1, j, k)], m_grids->omg_y[POS(i + 1, j, k)], m_grids->omg_z[POS(i + 1, j, k)]).norm();
        q = Vec3(m_grids->omg_x[POS(i - 1, j, k)], m_grids->omg_y[POS(i - 1, j, k)], m_grids->omg_z[POS(i - 1, j, k)]).norm();
        double grad1 = (p - q) * 0.5 / VOXEL_SIZE;

        p = Vec3(m_grids->omg_x[POS(i, j + 1, k)], m_grids->omg_y[POS(i, j + 1, k)], m_grids->omg_z[POS(i, j + 1, k)]).norm();
        q = Vec3(m_grids->omg_x[POS(i, j - 1, k)], m_grids->omg_y[POS(i, j - 1, k)], m_grids->omg_z[POS(i, j - 1, k)]).norm();
        double grad2 = (p - q) * 0.5 / VOXEL_SIZE;

        p = Vec3(m_grids->omg_x[POS(i, j, k + 1)], m_grids->omg_y[POS(i + 1, j, k + 1)], m_grids->omg_z[POS(i + 1, j, k + 1)]).norm();
        q = Vec3(m_grids->omg_x[POS(i, j, k - 1)], m_grids->omg_y[POS(i - 1, j, k - 1)], m_grids->omg_z[POS(i - 1, j, k - 1)]).norm();
        double grad3 = (p - q) * 0.5 / VOXEL_SIZE;

        Vec3 gradVort(grad1, grad2, grad3);
        // compute N vector
        Vec3 N_ijk = gradVort / (gradVort.norm() + 10e-20);

        Vec3 vorticity = Vec3(m_grids->omg_x[POS(i, j, k)], m_grids->omg_y[POS(i, j, k)], m_grids->omg_z[POS(i, j, k)]);
        Vec3 f = VORT_EPS * VOXEL_SIZE * vorticity.cross(N_ijk);
        m_grids->fx[POS(i, j, k)] += f[0];
        m_grids->fy[POS(i, j, k)] += f[1];
        m_grids->fz[POS(i, j, k)] += f[2];
    }
}

void Simulator::addForce()
{
    FOR_EACH_CELL
    {
        // ignore first cells
        if (i == 0 || j == 0 || k == 0)
        {
            continue;
        }

        m_grids->u[POSU(i, j, k)] += DT * (m_grids->fx[POS(i, j, k)] + m_grids->fx[POS(i + 1, j, k)]) * 0.5;
        m_grids->v[POSV(i, j, k)] += DT * (m_grids->fy[POS(i, j, k)] + m_grids->fy[POS(i, j + 1, k)]) * 0.5;
        m_grids->w[POSW(i, j, k)] += DT * (m_grids->fz[POS(i, j, k)] + m_grids->fz[POS(i, j, k + 1)]) * 0.5;

        m_grids->u0[POSU(i, j, k)] = m_grids->u[POSU(i, j, k)];
        m_grids->v0[POSV(i, j, k)] = m_grids->v[POSV(i, j, k)];
        m_grids->w0[POSW(i, j, k)] = m_grids->w[POSW(i, j, k)];
    }
}

void Simulator::advectVelocity()
{
    for (int k = 0; k <= N; ++k)
    {
        for (int j = 0; j <= N; ++j)
        {
            for (int i = 0; i <= N; ++i)
            {
                double x = i * LENGTH / (double)N;
                double y = (j + 0.5) * LENGTH / (double)N;
                double z = (k + 0.5) * LENGTH / (double)N;

                x = x - DT * macInterp(x, y, z, m_grids->u0, E_U, N + 1, N, N);
                y = y - DT * macInterp(x, y, z, m_grids->v0, E_V, N, N + 1, N);
                z = z - DT * macInterp(x, y, z, m_grids->w0, E_W, N, N, N + 1);

                m_grids->u[POSU(i, j, k)] = macInterp(x, y, z, m_grids->u0, E_U, N + 1, N, N);
            }
        }
    }
    for (int k = 0; k <= N; ++k)
    {
        for (int j = 0; j <= N; ++j)
        {
            for (int i = 0; i <= N; ++i)
            {
                double x = (i + 0.5) * LENGTH / (double)N;
                double y = j * LENGTH / (double)N;
                double z = (k + 0.5) * LENGTH / (double)N;

                x = x - DT * macInterp(x, y, z, m_grids->u0, E_U, N + 1, N, N);
                y = y - DT * macInterp(x, y, z, m_grids->v0, E_V, N, N + 1, N);
                z = z - DT * macInterp(x, y, z, m_grids->w0, E_W, N, N, N + 1);

                m_grids->v[POSV(i, j, k)] = macInterp(x, y, z, m_grids->v0, E_V, N, N + 1, N);
            }
        }
    }
    for (int k = 0; k <= N; ++k)
    {
        for (int j = 0; j <= N; ++j)
        {
            for (int i = 0; i <= N; ++i)
            {
                double x = (i + 0.5) * LENGTH / (double)N;
                double y = (j + 0.5) * LENGTH / (double)N;
                double z = k * LENGTH / (double)N;

                x = x - DT * macInterp(x, y, z, m_grids->u0, E_U, N + 1, N, N);
                y = y - DT * macInterp(x, y, z, m_grids->v0, E_V, N, N + 1, N);
                z = z - DT * macInterp(x, y, z, m_grids->w0, E_W, N, N, N + 1);

                m_grids->w[POSW(i, j, k)] = macInterp(x, y, z, m_grids->w0, E_W, N, N, N + 1);
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
                U[0] = (k > 0) ? m_grids->w[POSW(i, j, k - 1)] : 0.0;
                U[1] = (j > 0) ? m_grids->v[POSV(i, j - 1, k)] : 0.0;
                U[2] = (i > 0) ? m_grids->u[POSU(i - 1, j, k)] : 0.0;
                U[3] = m_grids->u[POSU(i + 1, j, k)];
                U[4] = m_grids->v[POSV(i, j + 1, k)];
                U[5] = m_grids->w[POSW(i, j, k + 1)];
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

    Eigen::Map<Eigen::VectorXd>(m_grids->pressure, SIZE) = x;
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
                    m_grids->u[POSU(i, j, k)] -= DT * (m_grids->pressure[POS(i + 1, j, k)] - m_grids->pressure[POS(i, j, k)]) * N / (double)LENGTH;
                }
                if (j < N - 1)
                {
                    m_grids->v[POSV(i, j, k)] -= DT * (m_grids->pressure[POS(i, j + 1, k)] - m_grids->pressure[POS(i, j, k)]) * N / (double)LENGTH;
                }
                if (k < N - 1)
                {
                    m_grids->w[POSW(i, j, k)] -= DT * (m_grids->pressure[POS(i, j, k + 1)] - m_grids->pressure[POS(i, j, k)]) * N / (double)LENGTH;
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

                x = x - DT * interp(x, y, z, m_grids->avg_u, N, N, N);
                y = y - DT * interp(x, y, z, m_grids->avg_v, N, N, N);
                z = z - DT * interp(x, y, z, m_grids->avg_w, N, N, N);

                m_grids->density[POS(i, j, k)] = interp(x, y, z, m_grids->density, N, N, N);
                m_grids->temperature[POS(i, j, k)] = interp(x, y, z, m_grids->temperature, N, N, N);
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

double Simulator::macInterp(double x, double y, double z, double q[], EMode mode, unsigned int Nx, unsigned int Ny, unsigned int Nz)
{
    x = std::fmax(0.0, std::fmin(Nx - 1 - 1e-6, N * x / (double)LENGTH));
    y = std::fmax(0.0, std::fmin(Ny - 1 - 1e-6, N * y / (double)LENGTH));
    z = std::fmax(0.0, std::fmin(Nz - 1 - 1e-6, N * z / (double)LENGTH));

    unsigned int i = x;
    unsigned int j = y;
    unsigned int k = z;
    double f[8];

    switch (mode)
    {
    case E_U:
        f[0] = q[POSU(i, j, k)];
        f[1] = q[POSU(i, j, k + 1)];
        f[2] = q[POSU(i, j + 1, k)];
        f[3] = q[POSU(i + 1, j, k)];
        f[4] = q[POSU(i, j + 1, k + 1)];
        f[5] = q[POSU(i + 1, j, k + 1)];
        f[6] = q[POSU(i + 1, j + 1, k)];
        f[7] = q[POSU(i + 1, j + 1, k + 1)];
        break;
    case E_V:
        f[0] = q[POSV(i, j, k)];
        f[1] = q[POSV(i, j, k + 1)];
        f[2] = q[POSV(i, j + 1, k)];
        f[3] = q[POSV(i + 1, j, k)];
        f[4] = q[POSV(i, j + 1, k + 1)];
        f[5] = q[POSV(i + 1, j, k + 1)];
        f[6] = q[POSV(i + 1, j + 1, k)];
        f[7] = q[POSV(i + 1, j + 1, k + 1)];
        break;
    case E_W:
        f[0] = q[POSW(i, j, k)];
        f[1] = q[POSW(i, j, k + 1)];
        f[2] = q[POSW(i, j + 1, k)];
        f[3] = q[POSW(i + 1, j, k)];
        f[4] = q[POSW(i, j + 1, k + 1)];
        f[5] = q[POSW(i + 1, j, k + 1)];
        f[6] = q[POSW(i + 1, j + 1, k)];
        f[7] = q[POSW(i + 1, j + 1, k + 1)];
        break;
    }

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
