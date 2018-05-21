#pragma once
#include "MACGrid.hpp"
#include "Volume.hpp"

class Scene
{
public:
  Scene(MACGrid *grids);
  ~Scene();

  void writeData();

  void update();
  void render();

private:
  void writeData_inVtiFormat();

  int m_file_num;
  MACGrid *m_grids;
  Volume *m_volume;
};
