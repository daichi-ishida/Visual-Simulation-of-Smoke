#pragma once
#include <string>
#include <memory>

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
  Volume(std::shared_ptr<MACGrid> grids);
  ~Volume();

  void initialize();

  void update();
  void draw() const;

  GLuint getProgramID() const;
  GLuint getVaoID() const;
  GLuint getCamPosID() const;
  GLuint getLightPosID() const;
  GLuint getLightIntensityID() const;
  GLuint getMatrixID() const;

private:
  void initVAO();
  void initShaders();

  std::string ReadFile(const std::string &filename);

  std::shared_ptr<MACGrid> m_grids;

  // ID
  GLuint vaoID;
  GLuint vboID;
  GLuint textureID;
  GLuint indexID;
  GLuint programID;
  GLuint vertexShaderID;
  GLuint fragmentShaderID;

  GLuint volumeTexID;

  GLuint cameraPosID;
  GLuint LightPosID;
  GLuint LightIntensityID;
  GLuint MatrixID;
  GLuint absorptionID;

  // data
  float absorption;
};