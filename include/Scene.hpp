#pragma once
#include <string>
#include <memory>

#define GLFW_INCLUDE_GLU
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

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
  Scene(std::shared_ptr<MACGrid> grids);
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

  std::shared_ptr<MACGrid> m_grids;
  std::unique_ptr<Camera> m_camera;
  std::unique_ptr<Volume> m_volume;
  std::unique_ptr<Wireframe> m_wireframe;

  std::unique_ptr<cv::VideoWriter> m_writer;
  std::string file_name;
};
