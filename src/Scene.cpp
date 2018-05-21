#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>

#define GLFW_INCLUDE_GLU
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include "Scene.hpp"
#include "constants.hpp"

Scene::Scene(MACGrid *grids) : m_file_num(0), m_grids(grids)
{
    // std::string vertex_shader_file = std::string("./shader/") + "volume.vert";
    // std::string fragment_shader_file = std::string("./shader/") + "volume.frag";
    // m_volume = new Volume(m_grids, vertex_shader_file, fragment_shader_file);
}

Scene::~Scene()
{
}

void Scene::writeData()
{
    writeData_inVtiFormat();
    ++m_file_num;
}

void Scene::update()
{
    m_volume->update();
}

void Scene::render()
{
    m_volume->draw();
}

/* private */
void Scene::writeData_inVtiFormat()
{
    std::ostringstream sout;
    sout << std::setfill('0') << std::setw(3) << std::right << m_file_num;
    std::string file_name;
    if (INTERPOLATION_METHOD == E_LINEAR)
    {
        file_name = "output/grids_linear" + sout.str() + ".vti";
    }
    else if (INTERPOLATION_METHOD == E_MONOTONIC_CUBIC)
    {
        file_name = "output/grids_monotonic" + sout.str() + ".vti";
    }

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
    ofs << "<ImageData WholeExtent='0 " << Nx << " 0 " << Ny << " 0 " << Nz << "' Origin='0 0 0' Spacing='" << VOXEL_SIZE << " " << VOXEL_SIZE << " " << VOXEL_SIZE << "'>" << std::endl;

    ofs
        << "<Piece Extent='0 " << Nx << " 0 " << Ny << " 0 " << Nz << "'>" << std::endl;

    ofs << "<CellData Vectors='velocity' Scalars='density temperature pressure vorticity'>" << std::endl;

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

    ofs << "<DataArray type='Float32' Name='vorticity' NumberOfComponents='1' format='ascii'>" << std::endl;
    for (int k = 0; k < Nz; ++k)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                ofs << m_grids->vort[POS(i, j, k)] << " ";
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