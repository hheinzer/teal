#include <math.h>

#include "simulation.h"
#include "teal/arena.h"
#include "teal/assert.h"
#include "teal/utils.h"

Simulation *simulation_create(const Equations *eqns, const char *prefix)
{
    assert(eqns);

    Simulation *sim = arena_calloc(1, sizeof(*sim));

    sim->eqns = eqns;
    sim->prefix = prefix;
    sim->courant = 0.99;  // NOLINT(readability-magic-numbers)

    sim->time.max = INFINITY;
    sim->time.output = INFINITY;

    sim->iter.max = NUMBER_MAX;
    sim->iter.output = NUMBER_MAX;

    RungeKutta ctx = {
        .time_order = eqns->space_order,
        .num_stages = eqns->space_order + 1,
    };
    simulation_set_advance(sim, "lserk", &ctx);

    return sim;
}
