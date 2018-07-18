#pragma once
#include <string>

#define GLFW_INCLUDE_GLU
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/noise.hpp>

#include "Camera.hpp"

class Wireframe
{
public:
  Wireframe();
  ~Wireframe();

  void initialize();

  void update();
  void draw() const;

  GLuint getProgramID() const;
  GLuint getVaoID() const;
  GLuint getMatrixID() const;

private:
  void initVAO();
  void initShaders();

  std::string ReadFile(const std::string &filename);

  // ID
  GLuint programID;
  GLuint vaoID;
  GLuint vboID;
  GLuint indexID;

  GLuint MatrixID;
};