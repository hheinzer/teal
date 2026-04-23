#include <assert.h>
#include <float.h>
#include <limits.h>
#include <string.h>

#include "simulation2.h"
#include "teal2.h"

Simulation *simulation2_create(const Equations *eqns, const char *prefix)
{
    assert(eqns && prefix);

    Simulation *sim = teal2_calloc(1, sizeof(*sim));
    sim->eqns = eqns;
    strcpy(sim->prefix, prefix);

    sim->time.max = DBL_MAX;
    sim->time.out = DBL_MAX;

    sim->iter.max = INT_MAX;
    sim->iter.out = INT_MAX;

    simulation2_set_advance(sim, (eqns->space_order == 1) ? "euler" : "lserk", 0.99, 0);

    return sim;
}
