#include <math.h>

#include "core/memory.h"
#include "core/sync.h"
#include "core/utils.h"
#include "equations.h"

void equations_residual(const Equations *eqns, double *residual)
{
    const long n_vars = eqns->n_vars;
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const double total_volume = eqns->mesh->volume;
    const ALIAS(cv, eqns->mesh->cell.volume);
    const double(*dudt)[n_vars] = (void *)eqns->vars.dudt;

    memory_setzero(residual, n_vars, sizeof(*residual));
    for (long i = 0; i < n_inner_cells; ++i)
        for (long v = 0; v < n_vars; ++v) residual[v] += cv[i] * sq(dudt[i][v]);
    for (long v = 0; v < n_vars; ++v) residual[v] = sqrt(sync_sum(residual[v]) / total_volume);
}
