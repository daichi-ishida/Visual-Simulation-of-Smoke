#include <iostream>

#define GLFW_INCLUDE_GLU
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include "constants.hpp"
#include "Scene.hpp"
#include "Simulator.hpp"
#include "MACGrid.hpp"

int main()
{
    if (glfwInit() == GL_FALSE)
    {
        fprintf(stderr, "Initialization failed!\n");
    }

    GLFWwindow *window = glfwCreateWindow(WIN_WIDTH, WIN_HEIGHT, WIN_TITLE,
                                          NULL, NULL);
    if (window == NULL)
    {
        fprintf(stderr, "Window creation failed!");
        glfwTerminate();
    }
    glfwMakeContextCurrent(window);
    if (glewInit() != GLEW_OK)
    {
        fprintf(stderr, "GLEW initialization failed!\n");
    }

    // set background color
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glEnable(GL_TEXTURE_3D);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);

    double time = 0.0;
    int step = 1;
    MACGrid *grids = new MACGrid();
    Simulator *simulator = new Simulator(grids, time);
    // simulator->update();
    Scene *scene = new Scene(grids);

    std::cout << "\n*** START SIMULATION ***\n";

    scene->writeData();
    scene->render();

    while (glfwWindowShouldClose(window) == GL_FALSE && glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS)
    // while (1)
    {
        std::cout << "\n=== STEP " << step << " ===\n";
        time += DT;
        simulator->update();

        scene->update();
        scene->writeData();
        scene->render();
        ++step;

        if (time >= FINISH_TIME)
        {
            break;
        }
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    std::cout << "\n*** END ***\n";

    if (simulator)
    {
        delete simulator;
    }
    if (scene)
    {
        delete scene;
    }
    if (grids)
    {
        delete grids;
    }

    return 0;
}