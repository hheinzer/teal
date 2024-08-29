#include <float.h>
#include <limits.h>
#include <string.h>

#include "restart.h"
#include "simulation.h"
#include "teal/memory.h"
#include "teal/option.h"

Simulation *simulation_create(Equations *eqns, const char *prefix)
{
    memory_sum_setzero();

    Simulation *sim = memory_calloc(1, sizeof(*sim));
    sim->eqns = eqns;
    strcpy(sim->prefix, prefix);

    simulation_set_max_time(sim, DBL_MAX);
    simulation_set_output_time(sim, DBL_MAX);
    simulation_set_max_iter(sim, LONG_MAX);
    simulation_set_output_iter(sim, LONG_MAX);
    simulation_set_abort(sim, -1, 0);
    simulation_set_cfl(sim, 0.99);
    simulation_set_advance(sim, "lserk");

    simulation_restart(sim, option.restart);

    return sim;
}
