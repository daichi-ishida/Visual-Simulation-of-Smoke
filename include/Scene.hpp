#pragma once
#include "Voxels.hpp"

class Scene
{
public:
  Scene(Voxels *voxels);
  ~Scene();

  void writeData();

private:
  void writeData_inVtiFormat();

  int m_file_num;
  Voxels *m_voxels;
};
