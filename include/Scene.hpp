#pragma once
#include <string>

#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/videoio.hpp>

#include "MACGrid.hpp"
#include "Volume.hpp"
#include "Camera.hpp"
#include "Wireframe.hpp"

class Scene
{
public:
  Scene(MACGrid *grids);
  ~Scene();

  void initialize();
  void update();
  void render();

  void writeData();

private:
  void writeData_inVtiFormat();
  void saveMovie();

  int m_file_num;
  glm::vec3 lightPos;
  glm::vec3 lightIntensity;

  MACGrid *m_grids;
  Camera *m_camera;
  Volume *m_volume;
  Wireframe *m_wireframe;

  cv::VideoWriter *m_writer;
  std::string file_name;
};
