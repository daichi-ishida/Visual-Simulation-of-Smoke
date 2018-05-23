#pragma once
#include <string>

#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/videoio.hpp>

#include "MACGrid.hpp"
#include "Volume.hpp"
#include "Camera.hpp"

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
  void saveMovie();

  int m_file_num;
  MACGrid *m_grids;
  Volume *m_volume;

  cv::VideoWriter *m_writer;
  std::string file_name;
};
