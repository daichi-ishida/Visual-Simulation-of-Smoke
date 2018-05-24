#pragma once
#define GLFW_INCLUDE_GLU
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/noise.hpp>

class Camera
{
public:
  Camera();
  ~Camera();

  void update();

  void GridViewControll();
  void FPScontroll();

  glm::vec3 getPos();
  glm::mat4 getProjectionMat();
  glm::mat4 getViewMat();

private:
  glm::vec3 m_position;
  glm::mat4 m_projectionMatix;
  glm::mat4 m_viewMatrix;
  glm::mat4 m_MVP;

  float m_lastTime;
  float m_currentTime;

  float m_horizontalAngle;
  float m_verticalAngle;
  float m_FoV;

  float m_speed;
  float m_mouseSpeed;
};