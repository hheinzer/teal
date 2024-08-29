#include "teal/print.h"

#include <float.h>
#include <limits.h>
#include <stdio.h>

#include "simulation.h"
#include "teal/memory.h"
#include "teal/option.h"
#include "teal/sync.h"
#include "teal/utils.h"

void simulation_print(const Simulation *sim)
{
    if (option.quiet) return;

    const alias(name, sim->eqns->vars.name);
    const char *abort = (sim->abort_variable < 0 ? "maximum" : name[sim->abort_variable]);
    const double dt = sim->cfl * equations_timestep(sim->eqns);

    const double size = sync_sum(memory_sum_get());

    if (sync.rank == 0) {
        printf("Simulation summary:\n");

        print_key("prefix", "%s", sim->prefix);

        if (sim->time > 0) print_key("initial time", "%g", sim->time);
        if (sim->max_time < DBL_MAX) print_key("max time", "%g", sim->max_time);
        if (sim->output_time < DBL_MAX) print_key("output time", "%g", sim->output_time);

        if (sim->max_iter < LONG_MAX) print_key("max iter", "%ld", sim->max_iter);
        if (sim->output_iter < LONG_MAX) print_key("output iter", "%ld", sim->output_iter);

        if (sim->abort_residual > 0) {
            print_key("abort variable", "%s", abort);
            print_key("abort residual", "%g", sim->abort_residual);
        }

        print_key("advance method", "%s", sim->advance.name);
        print_key("cfl number", "%g", sim->cfl);
        if (sim->time_order > 0) print_key("time order", "%ld", sim->time_order);
        if (sim->n_stages > 0) print_key("number of stages", "%ld", sim->n_stages);
        if (sim->tol_newton > 0) print_key("Newton tolerance", "%g", sim->tol_newton);
        if (sim->tol_krylov > 0) print_key("Krylov tolerance", "%g", sim->tol_krylov);
        if (sim->dim_krylov > 0) print_key("Krylov dimension", "%ld", sim->dim_krylov);

        print_key("initial dt", "%g", dt);

        print_size(size);
    }
}
