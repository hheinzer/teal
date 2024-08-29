#include <math.h>
#include <stdio.h>

#include "simulation.h"
#include "teal/option.h"
#include "teal/print.h"
#include "teal/sync.h"
#include "teal/utils.h"

double simulation_error(const Simulation *sim, long variable)
{
    const long n_inner_cells = sim->eqns->mesh->n_inner_cells;
    const double volume = sim->eqns->mesh->volume;
    const long n_vars = sim->eqns->n_vars;
    const long n_user = sim->eqns->n_user;
    const double time = sim->time;
    const alias(x, sim->eqns->mesh->cell.center);
    const alias(cv, sim->eqns->mesh->cell.volume);
    const double(*vars)[n_vars] = (void *)sim->eqns->vars.u;
    alias(compute, sim->eqns->user.compute);
    double user[n_user];

    double error = 0;
    for (long i = 0; i < n_inner_cells; ++i) {
        compute(user, vars[i], x[i], time);
        error += cv[i] * sq(user[variable] - vars[i][variable]);
    }
    error = sqrt(sync_sum(error) / volume);

    if (sync.rank == 0 && !option.quiet) {
        String key;
        sprintf(key, "L2 error %s", sim->eqns->vars.name[variable]);
        print_key(key, "%g", error);
    }

    return error;
}
