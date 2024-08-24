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
        if (sim->time_order > 0) printf(" | " KEYFMT ": %ld\n", "time order", sim->time_order);
        if (sim->n_stages > 0) printf(" | " KEYFMT ": %ld\n", "number of stages", sim->n_stages);
        if (sim->n_newton > 0)
            printf(" | " KEYFMT ": %ld\n", "number of Newton iterations", sim->n_newton);
        if (sim->n_krylov > 0)
            printf(" | " KEYFMT ": %ld\n", "number of Krylov iterations", sim->n_krylov);
        if (sim->tol_newton > 0) printf(" | " KEYFMT ": %g\n", "Newton tolerance", sim->tol_newton);
        if (sim->tol_krylov > 0) printf(" | " KEYFMT ": %g\n", "Krylov tolerance", sim->tol_krylov);
        printf(" | " KEYFMT ": %s\n", "advance function", sim->advance.name);
        printf(" | " KEYFMT ": %g\n", "cfl number", sim->cfl);
        if (sim->time > 0) printf(" | " KEYFMT ": %g\n", "initial time", sim->time);
        if (sim->iter > 0) printf(" | " KEYFMT ": %ld\n", "initial iter", sim->iter);
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

    const long n_cons = sim->eqns->n_cons;
    const long n_inner_cells = sim->eqns->mesh->n_inner_cells;
    if (sim->time_order != 1) size += n_inner_cells * n_cons * sizeof(*sim->advance.u0);

    const long n_krylov = sim->n_krylov;
    if (sim->time_order == 0) {
        size += n_inner_cells * n_cons * sizeof(sim->advance.xk);
        size += n_inner_cells * n_cons * sizeof(sim->advance.f0);
        size += n_inner_cells * n_cons * sizeof(sim->advance.fk);
        size += n_inner_cells * n_cons * sizeof(sim->advance.rk);
        size += n_inner_cells * n_cons * sizeof(sim->advance.dx);
        size += (n_krylov + 1) * n_inner_cells * n_cons * sizeof(sim->advance.V);
        size += n_krylov + 1 * sizeof(sim->advance.g);
        size += n_inner_cells * n_cons * sizeof(sim->advance.w);
        size += (n_krylov + 1) * n_krylov * sizeof(sim->advance.H);
        size += n_krylov * sizeof(sim->advance.s);
        size += n_krylov * sizeof(sim->advance.c);
        size += n_krylov * sizeof(sim->advance.y);
    }

    return size;
}
