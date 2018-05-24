
#define _USE_MATH_DEFINES
#include <cmath>
#include "Camera.hpp"
#include "constants.hpp"

Camera::Camera() : m_position(glm::vec3(0.0f, 0.0f, 5.0f)),
                   m_horizontalAngle(M_PI),
                   m_verticalAngle(0),
                   m_FoV(45.0f),
                   m_speed(7.0f),
                   m_mouseSpeed(0.0003f)
{
    update();
}

Camera::~Camera()
{
}

void Camera::update()
{
    // FPScontroll();
}

void Camera::GridViewControll()
{
    glm::vec3 center = glm::vec3((float)Nx * MAGNIFICATION / 2.0f, (float)Ny / 2.0f, (float)Nz / 2.0f);
    m_currentTime = glfwGetTime();
    float deltaTime = m_currentTime - m_lastTime;

    GLFWwindow *window = glfwGetCurrentContext();

    // Get mouse position
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);

    // Reset mouse position for next frame
    glfwSetCursorPos(window, WIN_WIDTH / 2, WIN_HEIGHT / 2); //reset mouse cursor position to the center wof window

    // Compute new orientation
    m_horizontalAngle += m_mouseSpeed * float(WIN_WIDTH / 2 - xpos);
    m_verticalAngle += m_mouseSpeed * float(WIN_HEIGHT / 2 - ypos);
    //restrict vertical angle not to be upside down
    if (m_verticalAngle < -M_PI / 2)
    {
        m_verticalAngle = -M_PI / 2;
    }
    if (m_verticalAngle > M_PI / 2)
    {
        m_verticalAngle = M_PI / 2;
    }

    // Direction : Spherical coordinates to Cartesian coordinates conversion
    glm::vec3 direction(
        std::cos(m_verticalAngle) * std::sin(m_horizontalAngle),
        std::sin(m_verticalAngle),
        std::cos(m_verticalAngle) * std::cos(m_horizontalAngle));

    // Right vector
    glm::vec3 right = glm::vec3(
        std::sin(m_horizontalAngle - M_PI / 2.0f), //substitute m_verticalAngle = 0 in the direction vector
        0,
        std::cos(m_horizontalAngle - M_PI / 2.0f));

    // Up vector
    glm::vec3 up = glm::cross(right, direction); //up vector which is perpendicular to right and direction , by definition

    // Move forward
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
    {
        m_position += direction * deltaTime * m_speed;
    }
    // Move backward
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
    {
        m_position -= direction * deltaTime * m_speed;
    }
    // Strafe right
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
    {
        m_position += right * deltaTime * m_speed;
    }
    // Strafe left
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
    {
        m_position -= right * deltaTime * m_speed;
    }
    //move up
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
    {
        m_position += up * deltaTime * m_speed;
    }
    //move down
    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
    {
        m_position -= up * deltaTime * m_speed;
    }

    // Projection matrix : 45ｰ Field of View,  ratio, display range : 0.1 unit <-> 100 units
    m_projectionMatix = glm::perspective(glm::radians(m_FoV), (float)WIN_WIDTH / (float)WIN_HEIGHT, 0.1f, 100.0f);

    // Camera matrix
    m_viewMatrix = glm::lookAt(
        m_position,             // Camera is here
        m_position + direction, // and looks here : at the same position, plus "direction"
        up                      // Head is up (set to 0,-1,0 to look upside-down)
    );

    // For the next frame, the "last time" will be "now"
    m_lastTime = glfwGetTime();
}

void Camera::FPScontroll()
{
    m_currentTime = glfwGetTime();
    float deltaTime = m_currentTime - m_lastTime;

    GLFWwindow *window = glfwGetCurrentContext();

    // Get mouse position
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);

    // Reset mouse position for next frame
    glfwSetCursorPos(window, WIN_WIDTH / 2, WIN_HEIGHT / 2); //reset mouse cursor position to the center wof window

    // Compute new orientation
    m_horizontalAngle += m_mouseSpeed * float(WIN_WIDTH / 2 - xpos);
    m_verticalAngle += m_mouseSpeed * float(WIN_HEIGHT / 2 - ypos);
    //restrict vertical angle not to be upside down
    if (m_verticalAngle < -M_PI / 2)
    {
        m_verticalAngle = -M_PI / 2;
    }
    if (m_verticalAngle > M_PI / 2)
    {
        m_verticalAngle = M_PI / 2;
    }

    // Direction : Spherical coordinates to Cartesian coordinates conversion
    glm::vec3 direction(
        std::cos(m_verticalAngle) * std::sin(m_horizontalAngle),
        std::sin(m_verticalAngle),
        std::cos(m_verticalAngle) * std::cos(m_horizontalAngle));

    // Right vector
    glm::vec3 right = glm::vec3(
        std::sin(m_horizontalAngle - M_PI / 2.0f), //substitute m_verticalAngle = 0 in the direction vector
        0,
        std::cos(m_horizontalAngle - M_PI / 2.0f));

    // Up vector
    glm::vec3 up = glm::cross(right, direction); //up vector which is perpendicular to right and direction , by definition

    // Move forward
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
    {
        m_position += direction * deltaTime * m_speed;
    }
    // Move backward
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
    {
        m_position -= direction * deltaTime * m_speed;
    }
    // Strafe right
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
    {
        m_position += right * deltaTime * m_speed;
    }
    // Strafe left
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
    {
        m_position -= right * deltaTime * m_speed;
    }
    //move up
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
    {
        m_position += up * deltaTime * m_speed;
    }
    //move down
    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
    {
        m_position -= up * deltaTime * m_speed;
    }

    // Projection matrix : 45ｰ Field of View,  ratio, display range : 0.1 unit <-> 100 units
    m_projectionMatix = glm::perspective(glm::radians(m_FoV), (float)WIN_WIDTH / (float)WIN_HEIGHT, 0.1f, 100.0f);

    // Camera matrix
    m_viewMatrix = glm::lookAt(
        m_position,             // Camera is here
        m_position + direction, // and looks here : at the same position, plus "direction"
        up                      // Head is up (set to 0,-1,0 to look upside-down)
    );

    // For the next frame, the "last time" will be "now"
    m_lastTime = glfwGetTime();
}

glm::vec3 Camera::getPos()
{
    return m_position;
}

glm::mat4 Camera::getProjectionMat()
{
    return m_projectionMatix;
}

glm::mat4 Camera::getViewMat()
{
    return m_viewMatrix;
}
