#include <cstdio>
#include "constants.hpp"

int main()
{
    double time = 0.0;
    printf("\n*** START SIMULATION ***\n");

    scene.writeData();

    while (1)
    {

        if (time >= FINISH_TIME)
        {
            break;
        }
    }

    printf("\n*** END ***\n");
    return 0;
}