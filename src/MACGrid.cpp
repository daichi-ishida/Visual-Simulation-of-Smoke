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
    double xstart = VOXEL_SIZE / 2.0;
    double ystart = VOXEL_SIZE / 2.0;
    double zstart = VOXEL_SIZE / 2.0;

    double x = xstart + i * VOXEL_SIZE;
    double y = ystart + j * VOXEL_SIZE;
    double z = zstart + k * VOXEL_SIZE;
    return Vec3(x, y, z);
}

Vec3 MACGrid::getVelocity(const Vec3 &pos)
{
    Vec3 vel;
    vel[0] = u0.interp(pos);
    vel[1] = v0.interp(pos);
    vel[2] = w0.interp(pos);
    return vel;
}
double MACGrid::getVelocityX(const Vec3 &pos)
{
    return u0.interp(pos);
}
double MACGrid::getVelocityY(const Vec3 &pos)
{
    return v0.interp(pos);
}
double MACGrid::getVelocityZ(const Vec3 &pos)
{
    return w0.interp(pos);
}

double MACGrid::getDensity(const Vec3 &pos)
{
    return density.interp(pos);
}

double MACGrid::getTemperature(const Vec3 &pos)
{
    return temperature.interp(pos);
}

double MACGrid::getPressure(const Vec3 &pos)
{
    return pressure.interp(pos);
}