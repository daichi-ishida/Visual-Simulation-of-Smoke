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
    Wireframe(Camera *camera);
    ~Wireframe();

    void update();
    void draw();

  private:
    std::string ReadFile(const std::string &filename);

    Camera *m_camera;

    // ID
    GLuint programID;
    GLuint vertexShaderID;
    GLuint fragmentShaderID;

    GLuint MatrixID;

    // data
    glm::mat4 MVP;

    GLuint m_vertex_array_object;
    GLuint m_vertex_buffer_object;
    GLuint m_index_buffer_object;
};