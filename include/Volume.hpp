#pragma once
#include <string>

#define GLFW_INCLUDE_GLU
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/noise.hpp>

#include "MACGrid.hpp"
#include "Camera.hpp"

class Volume
{
public:
  Volume(MACGrid *grids, const std::string &vertex_shader_file, const std::string &fragment_shader_file);
  ~Volume();

  void update();
  void draw();

private:
  std::string ReadFile(const std::string &filename);

  MACGrid *m_grids;
  Camera m_camera;

  // ID
  GLuint programID;
  GLuint vertexShaderID;
  GLuint fragmentShaderID;

  GLuint volumeTexID;
  GLuint cameraPosID;
  GLuint LightPosID;
  GLuint LightIntensityID;
  GLuint absorptionID;
  GLuint MatrixID;

  // data
  glm::vec3 cameraPos;
  glm::vec3 lightPos;
  glm::vec3 lightIntensity;
  float absorption;
  glm::mat4 MVP;

  GLuint m_vertex_array_object;
  GLuint m_vertex_buffer_object;
  GLuint m_texture_buffer_object;
  GLuint m_index_buffer_object;
};