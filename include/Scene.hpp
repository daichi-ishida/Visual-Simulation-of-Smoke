#pragma once
#include "MACGrid.hpp"

class Scene
{
public:
  Scene(MACGrid *grids);
  ~Scene();

  void writeData();

private:
  void writeData_inVtiFormat();

  int m_file_num;
  MACGrid *m_grids;
};
