#include <float.h>
#include <limits.h>
#include <stdio.h>

#include "core/sync.h"
#include "core/utils.h"
#include "simulation.h"
#include "teal.h"

static double compute_size(const Simulation *sim);

void simulation_print(const Simulation *sim)
{
    if (teal.quiet) return;

    const ALIAS(name, sim->eqns->vars.name);
    const char *var = (sim->abort_variable < 0 ? "maximum" : name[sim->abort_variable]);
    const double dt = sim->cfl * equations_timestep(sim->eqns);

    double size = sync_sum(compute_size(sim));
    const char mod = sizefmt(&size);

    if (teal.rank == 0) {
        printf("Simulation summary:\n");
        printf(" | " KEYFMT ": %s\n", "prefix", sim->prefix);
        if (teal.restart) printf(" | " KEYFMT ": %s\n", "restart", teal.restart);
        printf(" | " KEYFMT ": %ld\n", "time order", sim->time_order);
        printf(" | " KEYFMT ": %ld\n", "number of stages", sim->n_stages);
        printf(" | " KEYFMT ": %s\n", "advance function", sim->advance.name);
        printf(" | " KEYFMT ": %g\n", "cfl number", sim->cfl);
        if (sim->time > 0) printf(" | " KEYFMT ": %g\n", "time", sim->time);
        if (sim->iter > 0) printf(" | " KEYFMT ": %ld\n", "iter", sim->iter);
        if (sim->max_time < DBL_MAX) printf(" | " KEYFMT ": %g\n", "max time", sim->max_time);
        if (sim->max_iter < LONG_MAX) printf(" | " KEYFMT ": %ld\n", "max iter", sim->max_iter);
        if (sim->output_time < DBL_MAX)
            printf(" | " KEYFMT ": %g\n", "output time", sim->output_time);
        if (sim->output_iter < LONG_MAX)
            printf(" | " KEYFMT ": %ld\n", "output iter", sim->output_iter);
        if (sim->abort_residual > 0) {
            printf(" | " KEYFMT ": %s\n", "abort variable", var);
            printf(" | " KEYFMT ": %g\n", "abort residual", sim->abort_residual);
        }
        printf(" | " KEYFMT ": %g\n", "initial time step", dt);
        printf(" | " KEYFMT ": %g %cB\n", "memory size", size, mod);
    }
}

static double compute_size(const Simulation *sim)
{
    double size = sizeof(*sim);

    const long n_vars = sim->eqns->n_vars;
    const long n_inner_cells = sim->eqns->mesh->n_inner_cells;
    size += n_inner_cells * n_vars * sizeof(*sim->advance.buf);

    return size;
}
