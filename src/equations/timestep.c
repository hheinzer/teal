#include <float.h>

#include "core/sync.h"
#include "core/utils.h"
#include "equations.h"

double equations_timestep(const Equations *eqns)
{
    const long n_vars = eqns->n_vars;
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const ALIAS(cv, eqns->mesh->cell.volume);
    const ALIAS(p, eqns->mesh->cell.projection);
    ALIAS(time_step, eqns->timestep);
    const double(*u)[n_vars] = (void *)eqns->vars.u;

    double dt = DBL_MAX;
    for (long i = 0; i < n_inner_cells; ++i) dt = min(dt, time_step(eqns, u[i], p[i], cv[i]));
    return sync_min(dt);
}
