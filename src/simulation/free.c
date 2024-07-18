#include <stdlib.h>

#include "simulation.h"

void simulation_free(Simulation *sim)
{
    free(sim->advance.buf);

    *sim = (Simulation){0};
}
