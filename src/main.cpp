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
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
    glLineWidth(1.2f);

    glEnable(GL_TEXTURE_3D);
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE);
    // glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);

    double time = 0.0;
    int step = 1;
    MACGrid *grids = new MACGrid();
    Simulator *simulator = new Simulator(grids, time);
    Scene *scene = new Scene(grids);

    printf("\n*** START SIMULATION ***\n");

    // scene->writeData();
    scene->render();

    while (glfwWindowShouldClose(window) == GL_FALSE && glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS)
    // while (1)
    {
        printf("\n=== STEP %d ===\n", step);
        time += DT;
        simulator->update();

        scene->update();
        // scene->writeData();
        scene->render();
        ++step;

        if (time >= FINISH_TIME)
        {
            break;
        }
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    printf("\n*** END ***\n");

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