#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>

#define GLFW_INCLUDE_GLU
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include "Scene.hpp"
#include "constants.hpp"

Scene::Scene(MACGrid *grids) : m_file_num(0), m_grids(grids)
{
    m_camera = new Camera();
    m_volume = new Volume(m_grids, m_camera);
    m_wireframe = new Wireframe(m_camera);

    if (SAVE_MOVIE)
    {
        file_name = "output/grids_";
        switch (ADVECTION_METHOD)
        {
        case E_SEMI_LAGRANGE:
            file_name += "semi_lagrange_";
            break;
        case E_MAC_CORMACK:
            file_name += "mac_cormack_";
            break;
        }
        switch (INTERPOLATION_METHOD)
        {
        case E_LINEAR:
            file_name += "linear_";
            break;
        case E_MONOTONIC_CUBIC:
            file_name += "monotonic_";
            break;
        }
        file_name += ".avi";
        m_writer = new cv::VideoWriter(file_name, cv::VideoWriter::fourcc('M', 'J', 'P', 'G'), 30.0, cv::Size(WIN_WIDTH, WIN_HEIGHT));
        if (!m_writer->isOpened())
        {
            fprintf(stderr, "ERROR : file open error at writing data in .avi format\n %s cannot open\n", file_name.c_str());
            exit(EXIT_FAILURE);
        }
    }
}

Scene::~Scene()
{
    if (m_writer)
    {
        delete m_writer;
    }
    if (m_wireframe)
    {
        delete m_wireframe;
    }
    if (m_volume)
    {
        delete m_volume;
    }
    if (m_camera)
    {
        delete m_camera;
    }
}

void Scene::writeData()
{
    writeData_inVtiFormat();
    ++m_file_num;
}

void Scene::update()
{
    m_camera->update();
    m_volume->update();
    m_wireframe->update();
}

void Scene::render()
{
    m_volume->draw();
    m_wireframe->draw();
    if (SAVE_MOVIE)
    {
        saveMovie();
    }
}

/* private */
void Scene::writeData_inVtiFormat()
{
    std::ostringstream sout;
    sout << std::setfill('0') << std::setw(3) << std::right << m_file_num;
    file_name = "output/grids_";
    switch (ADVECTION_METHOD)
    {
    case E_SEMI_LAGRANGE:
        file_name += "semi_lagrange_";
        break;
    case E_MAC_CORMACK:
        file_name += "mac_cormack_";
        break;
    }
    switch (INTERPOLATION_METHOD)
    {
    case E_LINEAR:
        file_name += "linear_";
        break;
    case E_MONOTONIC_CUBIC:
        file_name += "monotonic_";
        break;
    }
    file_name += sout.str() + ".vti";

    std::ofstream ofs;
    ofs.open(file_name);
    if (!ofs)
    {
        fprintf(stderr, "ERROR : file open error at writing data in .vti format\n %s cannot open\n", file_name.c_str());
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

void Scene::saveMovie()
{
    glReadBuffer(GL_BACK);
    cv::Mat out_img(cv::Size(WIN_WIDTH, WIN_HEIGHT), CV_8UC3);
    glReadPixels(0, 0, WIN_WIDTH, WIN_HEIGHT, GL_BGR, GL_UNSIGNED_BYTE, out_img.data);
    cv::flip(out_img, out_img, 0);

    if (out_img.empty())
    {
        fprintf(stderr, "ERROR : no image %s \n", file_name.c_str());
        return;
    }

    *m_writer << out_img;
}