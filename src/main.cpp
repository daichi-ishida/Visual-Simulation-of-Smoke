#include <cstdio>
#include "constants.hpp"
#include "Scene.hpp"
#include "Simulator.hpp"
#include "Voxels.hpp"

int main()
{
    double time = 0.0;
    Voxels *voxels = new Voxels();
    Scene *scene = new Scene(voxels);
    Simulator *simulator = new Simulator(voxels);

    printf("\n*** START SIMULATION ***\n");

    //scene->writeData();

    while (1)
    {
        time += DT;
        simulator->update();
        //scene->writeData();

        if (time >= FINISH_TIME)
        {
            break;
        }
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
    if (voxels)
    {
        delete voxels;
    }

    return 0;
}