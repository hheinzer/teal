#include <math.h>

#include "equations.h"
#include "teal/memory.h"
#include "teal/sync.h"
#include "teal/utils.h"

void equations_residual(const Equations *eqns, double *residual)
{
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const double volume = eqns->mesh->volume;
    const long n_cons = eqns->n_cons;
    const alias(cv, eqns->mesh->cell.volume);
    const double(*dudt)[n_cons] = (void *)eqns->vars.dudt;

    memory_setzero(residual, n_cons, sizeof(*residual));
    for (long i = 0; i < n_inner_cells; ++i)
        for (long v = 0; v < n_cons; ++v) residual[v] += cv[i] * sq(dudt[i][v]);
    for (long v = 0; v < n_cons; ++v) residual[v] = sqrt(sync_sum(residual[v]) / volume);
}
