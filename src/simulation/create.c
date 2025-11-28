#include <assert.h>
#include <limits.h>
#include <math.h>

#include "simulation.h"
#include "teal/arena.h"

Simulation *simulation_create(const Equations *eqns, const char *prefix)
{
    assert(eqns);

    Simulation *sim = arena_calloc(1, sizeof(*sim));

    sim->eqns = eqns;
    sim->prefix = prefix;
    sim->courant = 0.99;  // NOLINT(readability-magic-numbers)

    sim->time.max = INFINITY;
    sim->time.out = INFINITY;

    sim->iter.max = LONG_MAX;
    sim->iter.out = LONG_MAX;

    simulation_set_advance(sim, (eqns->space_order == 1) ? "euler" : "lserk", 0);

    return sim;
}
