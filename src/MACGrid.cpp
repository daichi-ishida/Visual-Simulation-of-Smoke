#include "MACGrid.hpp"

MACGrid::MACGrid() : u(), v(), w(),
                     u0(), v0(), w0(),
                     avg_u(), avg_v(), avg_w(),
                     omg_x(), omg_y(), omg_z(), omg_length(),
                     eta_x(), eta_y(), eta_z(),
                     density(), temperature(), pressure(),
                     fx(), fy(), fz()
{
}

MACGrid::~MACGrid()
{
}