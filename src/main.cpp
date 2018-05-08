#include <iostream>
#include "constants.hpp"
#include "Scene.hpp"
#include "Simulator.hpp"
#include "MACGrid.hpp"

int main()
{
    double time = 0.0;
    int step = 1;
    MACGrid *grids = new MACGrid();
    Scene *scene = new Scene(grids);
    Simulator *simulator = new Simulator(grids, time);

    std::cout << "\n*** START SIMULATION ***\n";

    scene->writeData();

    while (1)
    {
        std::cout << "\n=== STEP " << step << " ===\n";
        time += DT;
        simulator->update();
        scene->writeData();
        ++step;

        if (time >= FINISH_TIME)
        {
            break;
        }
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