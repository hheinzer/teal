#include "equations.h"
#include "teal/memory.h"
#include "teal/sync.h"
#include "teal/utils.h"

void equations_average(const Equations *eqns, double *average)
{
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const double volume = eqns->mesh->volume;
    const long n_vars = eqns->n_vars;
    const alias(cv, eqns->mesh->cell.volume);
    const double(*u)[n_vars] = (void *)eqns->vars.u;

    memory_setzero(average, n_vars, sizeof(*average));
    for (long i = 0; i < n_inner_cells; ++i)
        for (long v = 0; v < n_vars; ++v) average[v] += u[i][v] * cv[i];
    for (long v = 0; v < n_vars; ++v) average[v] = sync_sum(average[v]) / volume;
}
