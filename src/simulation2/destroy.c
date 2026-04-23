#include <assert.h>

#include "simulation2.h"
#include "teal2.h"

void simulation2_destroy(Simulation *sim)
{
    assert(sim);

    teal2_free(sim->advance.context);

    teal2_free(sim);
}
