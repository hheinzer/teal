#include <float.h>

#include "equations.h"
#include "teal/sync.h"
#include "teal/utils.h"

double equations_timestep(Equations *eqns)
{
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const long n_vars = eqns->n_vars;
    const alias(cv, eqns->mesh->cell.volume);
    const alias(p, eqns->mesh->cell.projection);
    alias(compute, eqns->timestep.compute);
    const double(*u)[n_vars] = (void *)eqns->vars.u;
    double *dt = eqns->timestep.value;

    double min_dt = DBL_MAX;
    for (long i = 0; i < n_inner_cells; ++i) {
        dt[i] = compute(eqns, u[i], p[i], cv[i]);
        min_dt = min(min_dt, dt[i]);
    }
    return sync_min(min_dt);
}
