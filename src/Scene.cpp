#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

#include "constants.hpp"
#include "Scene.hpp"

Scene::Scene() : m_file_num(0)
{
}

Scene::~Scene()
{
}

void Scene::writeData()
{
}

/* private */

void Scene::writeData_inVtuFormat()
{
    std::ostringstream sout;
    sout << std::setfill('0') << std::setw(3) << std::right << m_file_num;

    std::string file_name = "output/particle_" + sout.str() + ".vtu";
    std::ofstream ofs;
    ofs.open(file_name);
    if (!ofs)
    {
        std::cout << "ERROR : file open error at writing data in .vtu format\n"
                  << file_name << " cannot open" << std::endl;
        exit(EXIT_FAILURE);
    }

    /* header */
    ofs << "<?xml version='1.0' encoding='UTF-8'?>" << std::endl;
    ofs << "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>" << std::endl;
    ofs << "<UnstructuredGrid>" << std::endl;

    ofs << "<Piece NumberOfCells='" << m_particle_set.getNumberOfParticles() << "' NumberOfPoints='" << m_particle_set.getNumberOfParticles() << "'>" << std::endl;

    /* point position */
    ofs << "<Points>" << std::endl;
    ofs << "<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>" << std::endl;
    for (unsigned int i = 0, n = m_particle_set.getNumberOfParticles(); i < n; ++i)
    {
        ofs << m_particle_set.positions[i].x() << " "
            << m_particle_set.positions[i].y() << " "
            << m_particle_set.positions[i].z() << std::endl;
    }
    ofs << "</DataArray>" << std::endl;
    ofs << "</Points>" << std::endl;

    /* point data */
    ofs << "<PointData>" << std::endl;

    // type
    ofs << "<DataArray NumberOfComponents='1' type='Int32' Name='ParticleType' format='ascii'>" << std::endl;
    for (unsigned int i = 0, n = m_particle_set.getNumberOfParticles(); i < n; ++i)
    {
        ofs << m_particle_set.types[i] << std::endl;
    }
    ofs << "</DataArray>" << std::endl;

    // velocity ( speed )
    ofs << "<DataArray NumberOfComponents='1' type='Float32' Name='Velocity' format='ascii'>" << std::endl;
    for (unsigned int i = 0, n = m_particle_set.getNumberOfParticles(); i < n; ++i)
    {
        double speed = m_particle_set.velocities[i].norm();
        ofs << speed << std::endl;
    }
    ofs << "</DataArray>" << std::endl;

    // pressures[i]
    ofs << "<DataArray NumberOfComponents='1' type='Float32' Name='Pressure' format='ascii'>" << std::endl;
    for (unsigned int i = 0, n = m_particle_set.getNumberOfParticles(); i < n; ++i)
    {
        ofs << m_particle_set.pressures[i] << std::endl;
    }
    ofs << "</DataArray>" << std::endl;

    // pressures[i] gradient
    ofs << "<DataArray NumberOfComponents='1' type='Float32' Name='Pressure Gradient' format='ascii'>" << std::endl;
    for (unsigned int i = 0, n = m_particle_set.getNumberOfParticles(); i < n; ++i)
    {
        double grad_magnitude = m_particle_set.pressure_gradients[i].norm();
        ofs << grad_magnitude << std::endl;
    }
    ofs << "</DataArray>" << std::endl;

    ofs << "</PointData>" << std::endl;

    ofs << "<Cells>" << std::endl;
    ofs << "<DataArray type='Int32' Name='connectivity' format='ascii'>" << std::endl;

    for (unsigned int i = 0, n = m_particle_set.getNumberOfParticles(); i < n; ++i)
    {
        ofs << i << std::endl;
    }
    ofs << "</DataArray>" << std::endl;
    ofs << "<DataArray type='Int32' Name='offsets' format='ascii'>" << std::endl;

    for (unsigned int i = 0, n = m_particle_set.getNumberOfParticles(); i < n; ++i)
    {
        ofs << i + 1 << std::endl;
    }
    ofs << "</DataArray>" << std::endl;
    ofs << "<DataArray type='UInt8' Name='types' format='ascii'>" << std::endl;

    for (unsigned int i = 0, n = m_particle_set.getNumberOfParticles(); i < n; ++i)
    {
        ofs << 1 << std::endl;
    }
    ofs << "</DataArray>" << std::endl;
    ofs << "</Cells>" << std::endl;

    ofs << "</Piece>" << std::endl;

    ofs << "</UnstructuredGrid>" << std::endl;
    ofs << "</VTKFile>" << std::endl;

    ofs.close();
}