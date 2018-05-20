#include <cmath>
#include <random>
#include <iostream>
#include "Simulator.hpp"

Simulator::Simulator(MACGrid *grids, double &time) : m_grids(grids), m_time(time), A(SIZE, SIZE), b(SIZE), x(SIZE)
{
    // nnz size is estimated by 7*SIZE because there are 7 nnz elements in a row.(center and neighbor 6)
    tripletList.reserve(7 * SIZE);
    ICCG.setTolerance(0.00001);

    /*set temperature */
    // std::random_device rnd;
    // std::mt19937 engine(rnd());
    // std::uniform_real_distribution<double> dist(0, T_AMP);
    OPENMP_FOR
    for (int k = 0; k < Nz; ++k)
    {
        OPENMP_FOR
        for (int j = 0; j < Ny; ++j)
        {
            OPENMP_FOR
            for (int i = 0; i < Nx; ++i)
            {
                m_grids->temperature(i, j, k) = (j / (float)Ny) * T_AMP + T_AMBIENT;
                // m_grids->temperature(i, j, k) = (j / (float)Ny) * T_AMP + dist(engine) + T_AMBIENT;
            }
        }
    }
    addSource();
    setEmitterVelocity();
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
    applyPressureTerm();
    advectScalar();
    if (m_time < EMIT_DURATION)
    {
        addSource();
        setEmitterVelocity();
    }
}

/* private */
void Simulator::addSource()
{
    OPENMP_FOR
    for (int k = Nz / 2 - SOURCE_SIZE_Z / 2; k < Nz / 2 + SOURCE_SIZE_Z / 2 + 1; ++k)
    {
        OPENMP_FOR
        for (int j = SOURCE_Y_MERGIN; j < SOURCE_Y_MERGIN + SOURCE_SIZE_Y; ++j)
        {
            OPENMP_FOR
            for (int i = Nx / 2 - SOURCE_SIZE_X / 2; i < Nx / 2 + SOURCE_SIZE_X / 2 + 1; ++i)
            {
                m_grids->density(i, j, k) = INIT_DENSITY;
            }
        }
    }
}

void Simulator::setEmitterVelocity()
{
    /* set emitter velocity */
    OPENMP_FOR
    for (int k = Nz / 2 - SOURCE_SIZE_Z / 2; k < Nz / 2 + SOURCE_SIZE_Z / 2 + 1; ++k)
    {
        OPENMP_FOR
        for (int j = SOURCE_Y_MERGIN; j < SOURCE_Y_MERGIN + SOURCE_SIZE_Y; ++j)
        {
            OPENMP_FOR
            for (int i = Nx / 2 - SOURCE_SIZE_X / 2; i < Nx / 2 + SOURCE_SIZE_X / 2 + 1; ++i)
            {
                m_grids->v(i, j, k) = INIT_VELOCITY;
                m_grids->v0(i, j, k) = m_grids->v(i, j, k);
            }
        }
    }
}

void Simulator::resetForce()
{

    OMP_FOR_EACH_CELL
    {
        m_grids->fx[POS(i, j, k)] = 0.0;
        m_grids->fy[POS(i, j, k)] = ALPHA * m_grids->density(i, j, k) - BETA * (m_grids->temperature(i, j, k) - T_AMBIENT);
        m_grids->fz[POS(i, j, k)] = 0.0;
    }
}

void Simulator::calVorticity()
{
    OMP_FOR_EACH_CELL
    {
        m_grids->avg_u[POS(i, j, k)] = (m_grids->u(i, j, k) + m_grids->u(i + 1, j, k)) * 0.5;
        m_grids->avg_v[POS(i, j, k)] = (m_grids->v(i, j, k) + m_grids->v(i, j + 1, k)) * 0.5;
        m_grids->avg_w[POS(i, j, k)] = (m_grids->w(i, j, k) + m_grids->w(i, j, k + 1)) * 0.5;
    }
    OMP_FOR_EACH_CELL
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

    OMP_FOR_EACH_CELL
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
        // compute gradient of vorticity
        double p, q;
        p = Vec3(m_grids->omg_x[POS(i + 1, j, k)], m_grids->omg_y[POS(i + 1, j, k)], m_grids->omg_z[POS(i + 1, j, k)]).norm();
        q = Vec3(m_grids->omg_x[POS(i - 1, j, k)], m_grids->omg_y[POS(i - 1, j, k)], m_grids->omg_z[POS(i - 1, j, k)]).norm();
        double grad1 = (p - q) * 0.5 / VOXEL_SIZE;

        p = Vec3(m_grids->omg_x[POS(i, j + 1, k)], m_grids->omg_y[POS(i, j + 1, k)], m_grids->omg_z[POS(i, j + 1, k)]).norm();
        q = Vec3(m_grids->omg_x[POS(i, j - 1, k)], m_grids->omg_y[POS(i, j - 1, k)], m_grids->omg_z[POS(i, j - 1, k)]).norm();
        double grad2 = (p - q) * 0.5 / VOXEL_SIZE;

        p = Vec3(m_grids->omg_x[POS(i, j, k + 1)], m_grids->omg_y[POS(i, j, k + 1)], m_grids->omg_z[POS(i, j, k + 1)]).norm();
        q = Vec3(m_grids->omg_x[POS(i, j, k - 1)], m_grids->omg_y[POS(i, j, k - 1)], m_grids->omg_z[POS(i, j, k - 1)]).norm();
        double grad3 = (p - q) * 0.5 / VOXEL_SIZE;

        Vec3 gradVort(grad1, grad2, grad3);
        // compute N vector
        Vec3 N_ijk(0, 0, 0);
        double norm = gradVort.norm();
        if (norm != 0)
        {
            N_ijk = gradVort / gradVort.norm();
        }

        Vec3 vorticity = Vec3(m_grids->omg_x[POS(i, j, k)], m_grids->omg_y[POS(i, j, k)], m_grids->omg_z[POS(i, j, k)]);
        Vec3 f = VORT_EPS * VOXEL_SIZE * vorticity.cross(N_ijk);
        m_grids->vort[POS(i, j, k)] = f.norm();
        m_grids->fx[POS(i, j, k)] += f[0];
        m_grids->fy[POS(i, j, k)] += f[1];
        m_grids->fz[POS(i, j, k)] += f[2];
    }
}

void Simulator::addForce()
{
    OMP_FOR_EACH_CELL
    {
        if (i < Nx - 1)
        {
            m_grids->u(i + 1, j, k) += DT * (m_grids->fx[POS(i, j, k)] + m_grids->fx[POS(i + 1, j, k)]) * 0.5;
        }
        if (j < Ny - 1)
        {
            m_grids->v(i, j + 1, k) += DT * (m_grids->fy[POS(i, j, k)] + m_grids->fy[POS(i, j + 1, k)]) * 0.5;
        }
        if (k < Nz - 1)
        {
            m_grids->w(i, j, k + 1) += DT * (m_grids->fz[POS(i, j, k)] + m_grids->fz[POS(i, j, k + 1)]) * 0.5;
        }

        m_grids->u0(i, j, k) = m_grids->u(i, j, k);
        m_grids->v0(i, j, k) = m_grids->v(i, j, k);
        m_grids->w0(i, j, k) = m_grids->w(i, j, k);
    }
}

void Simulator::advectVelocity()
{
    OMP_FOR_EACH_FACE_X
    {
        Vec3 pos_u = m_grids->getCenter(i, j, k) - 0.5 * Vec3(VOXEL_SIZE, 0, 0);
        Vec3 vel_u = m_grids->getVelocity(pos_u);
        Vec3 pos0_u = pos_u - DT * vel_u;
        m_grids->u(i, j, k) = m_grids->getVelocityX(pos0_u);
    }

    OMP_FOR_EACH_FACE_Y
    {
        Vec3 pos_v = m_grids->getCenter(i, j, k) - 0.5 * Vec3(0, VOXEL_SIZE, 0);
        Vec3 vel_v = m_grids->getVelocity(pos_v);
        Vec3 pos0_v = pos_v - DT * vel_v;
        m_grids->v(i, j, k) = m_grids->getVelocityY(pos0_v);
    }

    OMP_FOR_EACH_FACE_Z
    {
        Vec3 pos_w = m_grids->getCenter(i, j, k) - 0.5 * Vec3(0, 0, VOXEL_SIZE);
        Vec3 vel_w = m_grids->getVelocity(pos_w);
        Vec3 pos0_w = pos_w - DT * vel_w;
        m_grids->w(i, j, k) = m_grids->getVelocityZ(pos0_w);
    }
}

void Simulator::calPressure()
{
    tripletList.clear();
    A.setZero();
    b.setZero();
    x.setZero();

    double coeff = VOXEL_SIZE / DT;

#pragma omp parallel for ordered
    FOR_EACH_CELL
    {
        double F[6] = {k > 0, j > 0, i > 0, i < Nx - 1, j < Ny - 1, k < Nz - 1};
        double D[6] = {-1.0, -1.0, -1.0, 1.0, 1.0, 1.0};
        double U[6];
        U[0] = m_grids->w(i, j, k);
        U[1] = m_grids->v(i, j, k);
        U[2] = m_grids->u(i, j, k);
        U[3] = m_grids->u(i + 1, j, k);
        U[4] = m_grids->v(i, j + 1, k);
        U[5] = m_grids->w(i, j, k + 1);
        double sum_F = 0.0;

        for (int n = 0; n < 6; ++n)
        {
            sum_F += F[n];
            b(POS(i, j, k)) += D[n] * F[n] * U[n];
        }
        b(POS(i, j, k)) *= coeff;

#pragma omp ordered
        {
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

            if (i < Nx - 1)
            {
                tripletList.push_back(T(POS(i, j, k), POS(i + 1, j, k), F[3]));
            }
            if (j < Ny - 1)
            {
                tripletList.push_back(T(POS(i, j, k), POS(i, j + 1, k), F[4]));
            }
            if (k < Nz - 1)
            {
                tripletList.push_back(T(POS(i, j, k), POS(i, j, k + 1), F[5]));
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

    Eigen::Map<Eigen::VectorXd>(m_grids->pressure.getScalarPtr(), SIZE) = x;
}

void Simulator::applyPressureTerm()
{
    OMP_FOR_EACH_CELL
    {
        // compute gradient of pressure
        if (i < Nx - 1)
        {
            m_grids->u(i + 1, j, k) -= DT * (m_grids->pressure(i + 1, j, k) - m_grids->pressure(i, j, k)) / VOXEL_SIZE;
        }
        if (j < Ny - 1)
        {
            m_grids->v(i, j + 1, k) -= DT * (m_grids->pressure(i, j + 1, k) - m_grids->pressure(i, j, k)) / VOXEL_SIZE;
        }
        if (k < Nz - 1)
        {
            m_grids->w(i, j, k + 1) -= DT * (m_grids->pressure(i, j, k + 1) - m_grids->pressure(i, j, k)) / VOXEL_SIZE;
        }

        m_grids->u0(i, j, k) = m_grids->u(i, j, k);
        m_grids->v0(i, j, k) = m_grids->v(i, j, k);
        m_grids->w0(i, j, k) = m_grids->w(i, j, k);
    }
}

void Simulator::advectScalar()
{
    OMP_FOR_EACH_CELL
    {
        Vec3 pos_cell = m_grids->getCenter(i, j, k);
        Vec3 vel_cell = m_grids->getVelocity(pos_cell);
        Vec3 pos0_cell = pos_cell - DT * vel_cell;
        m_grids->density(i, j, k) = m_grids->getDensity(pos0_cell);
        m_grids->temperature(i, j, k) = m_grids->getTemperature(pos0_cell);
    }
}