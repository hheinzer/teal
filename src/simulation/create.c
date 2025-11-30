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

    sim->time.max = INFINITY;
    sim->time.out = INFINITY;

    sim->iter.max = LONG_MAX;
    sim->iter.out = LONG_MAX;

    const char *name = (eqns->space_order == 1) ? "euler" : "lserk";
    simulation_set_advance(sim, name, 0.99, 0);  // NOLINT(readability-magic-numbers)

    return sim;
}
