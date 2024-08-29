#include "free.h"

#include "advance.h"
#include "simulation.h"
#include "teal/memory.h"

void simulation_free(Simulation **sim)
{
    advance_free(*sim);
    memory_free(sim);
}

void advance_free(Simulation *sim)
{
    if (sim->advance.method == lserk)
        memory_free(&sim->advance.buf[0]);
    else if (sim->advance.method == implicit_euler)
        for (long i = 0; i < 12; ++i) memory_free(&sim->advance.buf[i]);
    memory_free(&sim->advance.buf);
}
