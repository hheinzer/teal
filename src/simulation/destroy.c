#include <assert.h>

#include "simulation.h"
#include "teal.h"

void simulation_destroy(Simulation *sim)
{
    assert(sim);

    teal_free(sim->advance.context);

    teal_free(sim);
}
