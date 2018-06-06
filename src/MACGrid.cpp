#include "MACGrid.hpp"

MACGrid::MACGrid() : avg_u(), avg_v(), avg_w(),
                     omg_x(), omg_y(), omg_z(), vort(),
                     fx(), fy(), fz()
{
}

MACGrid::~MACGrid()
{
}

Vec3 MACGrid::getCenter(int i, int j, int k)
{
    double half_dx = 0.5 * VOXEL_SIZE;

    double x = half_dx + i * VOXEL_SIZE;
    double y = half_dx + j * VOXEL_SIZE;
    double z = half_dx + k * VOXEL_SIZE;
    return Vec3(x, y, z);
}

Vec3 MACGrid::getVelocity(const Vec3 &pos)
{
    Vec3 vel;
    vel[0] = getVelocityX(pos);
    vel[1] = getVelocityY(pos);
    vel[2] = getVelocityZ(pos);
    return vel;
}
double MACGrid::getVelocityX(const Vec3 &pos)
{
    return u0.interp(pos - 0.5 * Vec3(0.0, VOXEL_SIZE, VOXEL_SIZE));
}
double MACGrid::getVelocityY(const Vec3 &pos)
{
    return v0.interp(pos - 0.5 * Vec3(VOXEL_SIZE, 0.0, VOXEL_SIZE));
}
double MACGrid::getVelocityZ(const Vec3 &pos)
{
    return w0.interp(pos - 0.5 * Vec3(VOXEL_SIZE, VOXEL_SIZE, 0.0));
}

double MACGrid::getDensity(const Vec3 &pos)
{
    return density0.interp(pos - 0.5 * Vec3(VOXEL_SIZE, VOXEL_SIZE, VOXEL_SIZE));
}

double MACGrid::getTemperature(const Vec3 &pos)
{
    return temperature0.interp(pos - 0.5 * Vec3(VOXEL_SIZE, VOXEL_SIZE, VOXEL_SIZE));
}

double MACGrid::getPressure(const Vec3 &pos)
{
    return pressure.interp(pos - 0.5 * Vec3(VOXEL_SIZE, VOXEL_SIZE, VOXEL_SIZE));
}