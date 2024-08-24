#include <math.h>

#include "core/memory.h"
#include "core/sync.h"
#include "core/utils.h"
#include "equations.h"

void equations_residual(const Equations *eqns, double *residual)
{
    const long n_cons = eqns->n_cons;
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const double total_volume = eqns->mesh->volume;
    const ALIAS(cv, eqns->mesh->cell.volume);
    const double(*dudt)[n_cons] = (void *)eqns->vars.dudt;

    memory_setzero(residual, n_cons, sizeof(*residual));
    for (long i = 0; i < n_inner_cells; ++i)
        for (long v = 0; v < n_cons; ++v) residual[v] += cv[i] * sq(dudt[i][v]);
    for (long v = 0; v < n_cons; ++v) residual[v] = sqrt(sync_sum(residual[v]) / total_volume);
}
