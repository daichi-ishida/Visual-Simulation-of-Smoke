#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

#include "constants.hpp"
#include "Scene.hpp"

Scene::Scene(MACGrid *grids) : m_file_num(0), m_grids(grids)
{
}

Scene::~Scene()
{
}

void Scene::writeData()
{
    writeData_inVtiFormat();
    ++m_file_num;
}

/* private */

void Scene::writeData_inVtiFormat()
{
    std::ostringstream sout;
    sout << std::setfill('0') << std::setw(3) << std::right << m_file_num;

    std::string file_name = "output/grids" + sout.str() + ".vti";
    std::ofstream ofs;
    ofs.open(file_name);
    if (!ofs)
    {
        std::cout << "ERROR : file open error at writing data in .vti format\n"
                  << file_name << " cannot open" << std::endl;
        exit(EXIT_FAILURE);
    }

    /* header */
    ofs << "<?xml version='1.0' encoding='UTF-8'?>" << std::endl;
    ofs << "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='ImageData'>" << std::endl;
    ofs << "<ImageData WholeExtent='0 " << Nx << " 0 " << Ny << " 0 " << Nz << "' Origin='0 0 0' Spacing='1.0 1.0 1.0'>>" << std::endl;

    ofs << "<Piece Extent='0 " << Nx << " 0 " << Ny << " 0 " << Nz << "'>" << std::endl;

    ofs << "<CellData Vectors='velocity' Scalars='density temperature pressure omega'>" << std::endl;

    ofs << "<DataArray type='Float32' Name='velocity' NumberOfComponents='3' format='ascii'>" << std::endl;
    FOR_EACH_CELL
    {
        ofs << m_grids->avg_u[POS(i, j, k)] << " " << m_grids->avg_v[POS(i, j, k)] << " " << m_grids->avg_w[POS(i, j, k)] << std::endl;
    }

    ofs << "</DataArray>" << std::endl;

    ofs << "<DataArray type='Float32' Name='density' NumberOfComponents='1' format='ascii'>" << std::endl;
    for (int k = 0; k < Nz; ++k)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                ofs << m_grids->density(i, j, k) << " ";
            }
            ofs << std::endl;
        }
    }
    ofs << "</DataArray>" << std::endl;

    ofs << "<DataArray type='Float32' Name='temperature' NumberOfComponents='1' format='ascii'>" << std::endl;
    for (int k = 0; k < Nz; ++k)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                ofs << m_grids->temperature(i, j, k) << " ";
            }
            ofs << std::endl;
        }
    }
    ofs << "</DataArray>" << std::endl;

    ofs << "<DataArray type='Float32' Name='pressure' NumberOfComponents='1' format='ascii'>" << std::endl;
    for (int k = 0; k < Nz; ++k)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                ofs << m_grids->pressure(i, j, k) << " ";
            }
            ofs << std::endl;
        }
    }
    ofs << "</DataArray>" << std::endl;

    ofs << "<DataArray type='Float32' Name='omega' NumberOfComponents='1' format='ascii'>" << std::endl;
    for (int k = 0; k < Nz; ++k)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                ofs << m_grids->omg_length[POS(i, j, k)] << " ";
            }
            ofs << std::endl;
        }
    }
    ofs << "</DataArray>" << std::endl;

    ofs << "</CellData>" << std::endl;
    ofs << "</Piece>" << std::endl;

    ofs << "</ImageData>" << std::endl;
    ofs << "</VTKFile>" << std::endl;

    ofs.close();
}