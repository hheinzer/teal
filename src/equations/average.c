#include "core/memory.h"
#include "core/sync.h"
#include "core/utils.h"
#include "equations.h"

void equations_average(const Equations *eqns, double *avg)
{
    const long n_vars = eqns->n_vars;
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const double volume = eqns->mesh->volume;
    const ALIAS(cv, eqns->mesh->cell.volume);
    const double(*u)[n_vars] = (void *)eqns->vars.u;

    memory_setzero(avg, n_vars, sizeof(*avg));
    for (long i = 0; i < n_inner_cells; ++i)
        for (long v = 0; v < n_vars; ++v) avg[v] += u[i][v] * cv[i];
    for (long v = 0; v < n_vars; ++v) avg[v] = sync_sum(avg[v]) / volume;
}
