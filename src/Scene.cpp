#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>

#include "Scene.hpp"
#include "constants.hpp"

Scene::Scene(std::shared_ptr<MACGrid> grids) : m_file_num(0), m_grids(grids)
{
    initialize();

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
        m_writer = std::make_unique<cv::VideoWriter>(file_name, cv::VideoWriter::fourcc('M', 'J', 'P', 'G'), 30.0, cv::Size(WIN_WIDTH, WIN_HEIGHT));
        if (!m_writer->isOpened())
        {
            fprintf(stderr, "ERROR : file open error at writing data in .avi format\n %s cannot open\n", file_name.c_str());
            exit(EXIT_FAILURE);
        }
    }
}

Scene::~Scene()
{
}

void Scene::initialize()
{
    m_camera = std::make_unique<Camera>();
    m_volume = std::make_unique<Volume>(m_grids);
    m_wireframe = std::make_unique<Wireframe>();
}

void Scene::update()
{
    float r = 7.0f;

    lightPos = glm::vec3(0.0f, -r, 0.0f);
    lightIntensity = 1.0f * glm::vec3(1.74f, 1.46f, 1.00f);
    // lightIntensity = glm::vec3(1.0f);

    m_camera->update();
    m_volume->update();
    m_wireframe->update();
}

void Scene::render()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glm::vec3 cameraPos = m_camera->getPos();
    GLuint programID = 0;

    /* --- draw volume --- */
    programID = m_volume->getProgramID();
    glUseProgram(programID);
    glBindVertexArray(m_volume->getVaoID());

    glUniform3f(m_volume->getLightPosID(), lightPos.x, lightPos.y, lightPos.z);
    glUniform3f(m_volume->getLightIntensityID(), lightIntensity.x, lightIntensity.y, lightIntensity.z);
    glUniform3f(m_volume->getCamPosID(), cameraPos.x, cameraPos.y, cameraPos.z);
    glUniformMatrix4fv(m_volume->getMatrixID(), 1, GL_FALSE, &(m_camera->getMVP())[0][0]);

    m_volume->draw();

    glBindVertexArray(0);
    glUseProgram(0);
    /* --- end --- */

    /* --- draw wireframe --- */
    programID = m_wireframe->getProgramID();
    glUseProgram(programID);
    glBindVertexArray(m_wireframe->getVaoID());

    glUniformMatrix4fv(m_wireframe->getMatrixID(), 1, GL_FALSE, &(m_camera->getMVP())[0][0]);

    m_wireframe->draw();

    glBindVertexArray(0);
    glUseProgram(0);
    /* --- end --- */

    if (SAVE_MOVIE)
    {
        saveMovie();
    }
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