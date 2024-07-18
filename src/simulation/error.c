#include <math.h>
#include <stdio.h>

#include "core/sync.h"
#include "core/utils.h"
#include "simulation.h"
#include "teal.h"

void simulation_error(const Simulation *sim, long ivars, long iuser)
{
    const double time = sim->time;
    const long n_vars = sim->eqns->n_vars;
    const long n_user = sim->eqns->n_user;
    const long n_inner_cells = sim->eqns->mesh->n_inner_cells;
    const double total_volume = sim->eqns->mesh->volume;
    const ALIAS(x, sim->eqns->mesh->cell.center);
    const ALIAS(v, sim->eqns->mesh->cell.volume);
    const double(*vars)[n_vars] = (void *)sim->eqns->vars.u;
    ALIAS(compute, sim->eqns->user.compute);
    double user[n_user];
    double error = 0;

    for (long i = 0; i < n_inner_cells; ++i) {
        compute(user, vars[i], x[i], time);
        error += v[i] * (user[iuser] - vars[i][ivars]) * (user[iuser] - vars[i][ivars]);
    }
    error = sqrt(sync_sum(error) / total_volume);

    if (teal.rank == 0) {
        char key[128];
        snprintf(key, sizeof(key), "L2 error %s", sim->eqns->vars.name[ivars]);
        printf(" | " KEYFMT ": %g\n", key, error);
    }
}
