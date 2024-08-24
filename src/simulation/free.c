#include <stdlib.h>

#include "simulation.h"

void simulation_free(Simulation *sim)
{
    free(sim->advance.u0);

    free(sim->advance.xk);
    free(sim->advance.f0);
    free(sim->advance.fk);
    free(sim->advance.rk);
    free(sim->advance.dx);
    free(sim->advance.V);
    free(sim->advance.g);
    free(sim->advance.w);
    free(sim->advance.H);
    free(sim->advance.s);
    free(sim->advance.c);
    free(sim->advance.y);

    *sim = (Simulation){0};
}
