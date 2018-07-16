#include <cmath>
#include <random>
#include "Simulator.hpp"

Simulator::Simulator(std::shared_ptr<MACGrid> grids, double &time) : m_grids(grids), m_time(time), A(SIZE, SIZE), b(SIZE), x(SIZE)
{
    // nnz size is estimated by 7*SIZE because there are 7 nnz elements in a row.(center and neighbor 6)
    tripletList.reserve(7 * SIZE);
    ICCG.setTolerance(1e-8);

    /*set temperature */
    std::random_device rnd;
    std::mt19937 engine(rnd());
    std::uniform_real_distribution<double> dist(0, T_AMP);

    OPENMP_FOR_COLLAPSE
    FOR_EACH_CELL
    {
        // m_grids->temperature(i, j, k) = (j / (float)Ny) * T_AMP + T_AMBIENT;
        m_grids->temperature(i, j, k) = (j / (float)Ny) * T_AMP + dist(engine) + T_AMBIENT;
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
    calPressure();
    applyPressureTerm();
    advectVelocity();
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
    switch (EMITTER_POS)
    {
    case E_TOP:
    {
        OPENMP_FOR_COLLAPSE
        for (int k = (Nz - SOURCE_SIZE_Z) / 2; k < (Nz + SOURCE_SIZE_Z) / 2; ++k)
        {
            for (int j = SOURCE_Y_MERGIN; j < SOURCE_Y_MERGIN + SOURCE_SIZE_Y; ++j)
            {
                for (int i = (Nx - SOURCE_SIZE_X) / 2; i < (Nx + SOURCE_SIZE_X) / 2; ++i)
                {
                    m_grids->density(i, j, k) = INIT_DENSITY;
                }
            }
        }
        break;
    }

    case E_BOTTOM:
    {
        OPENMP_FOR_COLLAPSE
        for (int k = (Nz - SOURCE_SIZE_Z) / 2; k < (Nz + SOURCE_SIZE_Z) / 2; ++k)
        {
            for (int j = Ny - SOURCE_Y_MERGIN - SOURCE_SIZE_Y; j < Ny - SOURCE_Y_MERGIN; ++j)
            {
                for (int i = (Nx - SOURCE_SIZE_X) / 2; i < (Nx + SOURCE_SIZE_X) / 2; ++i)
                {
                    m_grids->density(i, j, k) = INIT_DENSITY;
                }
            }
        }
        break;
    }
    }
}

void Simulator::setEmitterVelocity()
{
    switch (EMITTER_POS)
    {
    case E_TOP:
    {
        OPENMP_FOR_COLLAPSE
        for (int k = (Nz - SOURCE_SIZE_Z) / 2; k < (Nz + SOURCE_SIZE_Z) / 2; ++k)
        {
            for (int j = SOURCE_Y_MERGIN; j < SOURCE_Y_MERGIN + SOURCE_SIZE_Y; ++j)
            {
                for (int i = (Nx - SOURCE_SIZE_X) / 2; i < (Nx + SOURCE_SIZE_X) / 2; ++i)
                {
                    m_grids->v(i, j, k) = INIT_VELOCITY;
                    m_grids->v0(i, j, k) = m_grids->v(i, j, k);
                }
            }
        }
        break;
    }

    case E_BOTTOM:
    {
        OPENMP_FOR_COLLAPSE
        for (int k = (Nz - SOURCE_SIZE_Z) / 2; k < (Nz + SOURCE_SIZE_Z) / 2; ++k)
        {
            for (int j = Ny - SOURCE_Y_MERGIN - SOURCE_SIZE_Y; j < Ny - SOURCE_Y_MERGIN + 1; ++j)
            {
                for (int i = (Nx - SOURCE_SIZE_X) / 2; i < (Nx + SOURCE_SIZE_X) / 2; ++i)
                {
                    m_grids->v(i, j, k) = -INIT_VELOCITY;
                    m_grids->v0(i, j, k) = m_grids->v(i, j, k);
                }
            }
        }
        break;
    }
    }
}

void Simulator::resetForce()
{
    OPENMP_FOR_COLLAPSE
    FOR_EACH_CELL
    {
        m_grids->fx[POS(i, j, k)] = 0.0;
        m_grids->fy[POS(i, j, k)] = ALPHA * m_grids->density(i, j, k) - BETA * (m_grids->temperature(i, j, k) - T_AMBIENT);
        m_grids->fz[POS(i, j, k)] = 0.0;
    }
}

void Simulator::calVorticity()
{
    OPENMP_FOR_COLLAPSE
    FOR_EACH_CELL
    {
        m_grids->avg_u[POS(i, j, k)] = (m_grids->u(i, j, k) + m_grids->u(i + 1, j, k)) * 0.5;
        m_grids->avg_v[POS(i, j, k)] = (m_grids->v(i, j, k) + m_grids->v(i, j + 1, k)) * 0.5;
        m_grids->avg_w[POS(i, j, k)] = (m_grids->w(i, j, k) + m_grids->w(i, j, k + 1)) * 0.5;
    }

    OPENMP_FOR_COLLAPSE
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

    OPENMP_FOR_COLLAPSE
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
    OPENMP_FOR_COLLAPSE
    FOR_EACH_CELL
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
    }
}

void Simulator::calPressure()
{
    tripletList.clear();
    A.setZero();
    b.setZero();
    x.setZero();

    double coeff = VOXEL_SIZE / DT;

#pragma omp parallel for collapse(3) ordered
    FOR_EACH_CELL
    {
        double F[6] = {static_cast<double>(k > 0), static_cast<double>(j > 0), static_cast<double>(i > 0),
                       static_cast<double>(i < Nx - 1), static_cast<double>(j < Ny - 1), static_cast<double>(k < Nz - 1)};
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
        printf("SUCCESS: Convergence\n");
    }
    else
    {
        fprintf(stderr, "FAILED: No Convergence\n");
    }
    x = ICCG.solve(b);
    printf("#iterations:     %d \n", static_cast<int>(ICCG.iterations()));
    printf("estimated error: %e \n", ICCG.error());

    Eigen::Map<Eigen::VectorXd>(m_grids->pressure.begin(), SIZE) = x;
}

void Simulator::applyPressureTerm()
{
    OPENMP_FOR_COLLAPSE
    FOR_EACH_CELL
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
    }
    std::copy(m_grids->u.begin(), m_grids->u.end(), m_grids->u0.begin());
    std::copy(m_grids->v.begin(), m_grids->v.end(), m_grids->v0.begin());
    std::copy(m_grids->w.begin(), m_grids->w.end(), m_grids->w0.begin());
}

void Simulator::advectVelocity()
{
    switch (ADVECTION_METHOD)
    {
    case E_SEMI_LAGRANGE:
    {
        OPENMP_FOR_COLLAPSE
        FOR_EACH_FACE_X
        {
            Vec3 pos_u = m_grids->getCenter(i, j, k) - 0.5 * Vec3(VOXEL_SIZE, 0, 0);
            Vec3 vel_u = m_grids->getVelocity(pos_u);
            pos_u -= DT * vel_u;
            m_grids->u(i, j, k) = m_grids->getVelocityX(pos_u);
        }

        OPENMP_FOR_COLLAPSE
        FOR_EACH_FACE_Y
        {
            Vec3 pos_v = m_grids->getCenter(i, j, k) - 0.5 * Vec3(0, VOXEL_SIZE, 0);
            Vec3 vel_v = m_grids->getVelocity(pos_v);
            pos_v -= DT * vel_v;
            m_grids->v(i, j, k) = m_grids->getVelocityY(pos_v);
        }

        OPENMP_FOR_COLLAPSE
        FOR_EACH_FACE_Z
        {
            Vec3 pos_w = m_grids->getCenter(i, j, k) - 0.5 * Vec3(0, 0, VOXEL_SIZE);
            Vec3 vel_w = m_grids->getVelocity(pos_w);
            pos_w -= DT * vel_w;
            m_grids->w(i, j, k) = m_grids->getVelocityZ(pos_w);
        }
        break;
    }

    case E_MAC_CORMACK:
    {
        OPENMP_FOR_COLLAPSE
        FOR_EACH_FACE_X
        {
            double u_n = m_grids->u0(i, j, k);
            Vec3 pos_u = m_grids->getCenter(i, j, k) - 0.5 * Vec3(VOXEL_SIZE, 0, 0);
            Vec3 vel_u = m_grids->getVelocity(pos_u);
            // forward advection
            pos_u -= DT * vel_u;
            double u_np1_hat = m_grids->getVelocityX(pos_u);
            // backward advection
            pos_u += DT * m_grids->getVelocity(pos_u);
            double u_n_hat = m_grids->getVelocityX(pos_u);

            m_grids->u(i, j, k) = u_np1_hat + 0.5 * (u_n - u_n_hat);
        }
        OPENMP_FOR_COLLAPSE
        FOR_EACH_FACE_Y
        {
            double v_n = m_grids->v0(i, j, k);
            Vec3 pos_v = m_grids->getCenter(i, j, k) - 0.5 * Vec3(0, VOXEL_SIZE, 0);
            Vec3 vel_v = m_grids->getVelocity(pos_v);
            // forward advection
            pos_v -= DT * vel_v;
            double v_np1_hat = m_grids->getVelocityY(pos_v);
            // backward advection
            pos_v += DT * m_grids->getVelocity(pos_v);
            double v_n_hat = m_grids->getVelocityY(pos_v);

            m_grids->v(i, j, k) = v_np1_hat + 0.5 * (v_n - v_n_hat);
        }

        OPENMP_FOR_COLLAPSE
        FOR_EACH_FACE_Z
        {
            double w_n = m_grids->w0(i, j, k);
            Vec3 pos_w = m_grids->getCenter(i, j, k) - 0.5 * Vec3(0, 0, VOXEL_SIZE);
            Vec3 vel_w = m_grids->getVelocity(pos_w);
            // forward advection
            pos_w -= DT * vel_w;
            double w_np1_hat = m_grids->getVelocityZ(pos_w);
            // backward advection
            pos_w += DT * m_grids->getVelocity(pos_w);
            double w_n_hat = m_grids->getVelocityZ(pos_w);

            m_grids->w(i, j, k) = w_np1_hat + 0.5 * (w_n - w_n_hat);
        }
        break;
    }
    }
}

void Simulator::advectScalar()
{
    std::copy(m_grids->density.begin(), m_grids->density.end(), m_grids->density0.begin());
    std::copy(m_grids->temperature.begin(), m_grids->temperature.end(), m_grids->temperature0.begin());

    switch (ADVECTION_METHOD)
    {
    case E_SEMI_LAGRANGE:
    {
        OPENMP_FOR_COLLAPSE
        FOR_EACH_CELL
        {
            Vec3 pos_cell = m_grids->getCenter(i, j, k);
            Vec3 vel_cell = m_grids->getVelocity(pos_cell);
            pos_cell -= DT * vel_cell;
            m_grids->density(i, j, k) = m_grids->getDensity(pos_cell);
            m_grids->temperature(i, j, k) = m_grids->getTemperature(pos_cell);
        }
        break;
    }
    case E_MAC_CORMACK:
    {
        OPENMP_FOR_COLLAPSE
        FOR_EACH_CELL
        {
            double d_n = m_grids->density(i, j, k);
            double t_n = m_grids->temperature(i, j, k);
            Vec3 pos_cell = m_grids->getCenter(i, j, k);
            Vec3 vel_cell = m_grids->getVelocity(pos_cell);
            // forward advection
            pos_cell -= DT * vel_cell;
            double d_np1_hat = m_grids->getDensity(pos_cell);
            double t_np1_hat = m_grids->getTemperature(pos_cell);
            // backward advection
            pos_cell += DT * m_grids->getVelocity(pos_cell);
            double d_n_hat = m_grids->getDensity(pos_cell);
            double t_n_hat = m_grids->getTemperature(pos_cell);

            m_grids->density(i, j, k) = d_np1_hat + 0.5 * (d_n - d_n_hat);
            m_grids->temperature(i, j, k) = t_np1_hat + 0.5 * (t_n - t_n_hat);
        }
        break;
    }
    }
}