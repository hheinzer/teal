#include "simulation.h"

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "advance.h"
#include "array.h"
#include "equations.h"
#include "global.h"
#include "memory.h"
#include "mesh.h"
#include "sync.h"
#include "utils.h"

static char m_advance[128];  // advance method
static volatile sig_atomic_t m_terminate = 0;

static void terminate(int signum);
static void write(Simulation *sim, double *output_time, long *output_iter);

Simulation simulation_create(const char *prefix, Equations *eqns)
{
    Simulation sim = {
        .max_time = DBL_MAX,
        .output_time = DBL_MAX,
        .max_iter = LONG_MAX,
        .output_iter = LONG_MAX,
        .cfl = 0.99,
        .abort_variable = -1,
        .abort_residual = 0,
        .prefix = utils_strdup(prefix),
        .eqns = eqns,
    };
    char fname[128];
    snprintf(fname, sizeof(fname), "%s_residuals.dat", prefix);
    if (eqns->mesh->rank == 0) {
        const long n_vars = eqns->vars.n_fields;
        const ALIAS(name, eqns->vars.name);
        sim.residual_file = fopen(fname, "w");
        for (long v = 0; v < n_vars; ++v) fprintf(sim.residual_file, "%s ", name[v]);
        fputc('\n', sim.residual_file);
    }
    simulation_set_time_order(&sim, eqns->space_order, eqns->space_order + 1);
    return sim;
}

void simulation_free(Simulation *sim)
{
    if (sim->eqns->mesh->rank == 0) fclose(sim->residual_file);
    free(sim->buf);
    free(sim->prefix);
    *sim = (typeof(*sim)){};
}

void simulation_set_time_order(Simulation *sim, const long time_order, const long n_stages)
{
    assert(1 <= time_order && time_order <= 3 && "invalid time order");
    assert(1 <= n_stages && n_stages <= 6 && "invalid number of stages");
    sim->time_order = time_order;
    sim->n_stages = n_stages;
    if (time_order == 1 && n_stages == 1) {
        sim->advance = advance_euler;
        strncpy(m_advance, "euler", sizeof(m_advance) - 1);

        free(sim->buf);
        sim->buf = 0;
    }
    else {
        sim->advance = advance_lserk;
        strncpy(m_advance, "lserk", sizeof(m_advance) - 1);

        const long n_vars = sim->eqns->vars.n_fields;
        const long n_inner_cells = sim->eqns->mesh->n_inner_cells;
        sim->buf = memory_calloc(n_inner_cells * n_vars, sizeof(*sim->buf));
    }
}

void simulation_print(const Simulation *sim)
{
    const double dt = sync_min(sim->eqns->time_step(sim->eqns));

    if (sim->eqns->mesh->rank == 0) {
        printf("Simulation summary:\n");

        printf(" | " FMT_KEY ": %ld\n", "time order", sim->time_order);
        printf(" | " FMT_KEY ": %ld\n", "number of stages", sim->n_stages);
        printf(" | " FMT_KEY ": %s\n", "advance method", m_advance);

        if (sim->max_time < DBL_MAX) printf(" | " FMT_KEY ": %g\n", "max time", sim->max_time);
        if (sim->max_iter < LONG_MAX) printf(" | " FMT_KEY ": %ld\n", "max iter", sim->max_iter);

        if (sim->output_time < DBL_MAX)
            printf(" | " FMT_KEY ": %g\n", "output time", sim->output_time);
        if (sim->output_iter < LONG_MAX)
            printf(" | " FMT_KEY ": %ld\n", "output iter", sim->output_iter);

        printf(" | " FMT_KEY ": %g\n", "cfl", sim->cfl);

        if (sim->abort_residual > 0) {
            printf(
                " | " FMT_KEY ": %s\n", "abort variable",
                (sim->abort_variable < 0 ? "maximum" : sim->eqns->vars.name[sim->abort_variable]));
            printf(" | " FMT_KEY ": %g\n", "abort residual", sim->abort_residual);
        }

        printf(" | " FMT_KEY ": %g\n", "initial time step", sim->cfl * dt);
    }
}

void simulation_run(Simulation *sim)
{
    double output_time = sim->output_time;
    long output_iter = sim->output_iter;
    mesh_write(sim->eqns->mesh, sim->prefix);
    write(sim, &output_time, &output_iter);

    const long n_vars = sim->eqns->vars.n_fields;
    double residual[n_vars] = {};
    bool converged = false;

    signal(SIGINT, terminate);
    signal(SIGTERM, terminate);

    const bool rank = sim->eqns->mesh->rank;
    if (rank == 0)
        printf(" |  %13s  %13s  %13s  %13s  %13s\n", "iter", "time", "dt", "residual", "wtime");

    const double timer_start = MPI_Wtime();
    double timer_last = timer_start;

    while (sim->time < sim->max_time && sim->iter < sim->max_iter && !converged && !m_terminate) {
        const double max_dt = sim->cfl * sync_min(sim->eqns->time_step(sim->eqns));
        const double dt = MIN(max_dt, MIN(sim->max_time, output_time) - sim->time);
        assert(isfinite(dt) && "this can't be good");

        sim->advance(sim, dt);
        sim->time += dt;
        sim->iter += 1;

        equations_residual(sim->eqns, residual);
        if (rank == 0) array_fprint(sim->residual_file, residual, n_vars);
        const double max_residual =
            (sim->abort_variable < 0 ? array_max(residual, n_vars) : residual[sim->abort_variable]);
        converged = (max_residual < sim->abort_residual);

        if (sim->time >= MIN(sim->max_time, output_time) ||
            sim->iter >= MIN(sim->max_iter, output_iter) || converged || m_terminate) {
            write(sim, &output_time, &output_iter);
            const double timer_curr = MPI_Wtime();
            if (rank == 0)
                printf(" |  %13ld  %13g  %13g  %13g  %13g\n", sim->iter, sim->time, max_dt,
                       max_residual, timer_curr - timer_last);
            timer_last = timer_curr;
        }
    }

    const double timer_stop = MPI_Wtime();

    if (rank == 0) printf(" | " FMT_KEY ": %g\n", "computation time", timer_stop - timer_start);
}

void simulation_error(const Simulation *sim, const long i_vars, const long i_user)
{
    const long n_user = sim->eqns->user.n_fields;
    const long n_inner_cells = sim->eqns->mesh->n_inner_cells;
    const double total_volume = sim->eqns->mesh->volume;
    const ALIAS(x, sim->eqns->mesh->cell.center);
    const ALIAS(volume, sim->eqns->mesh->cell.volume);
    const FIELDS(vars, sim->eqns->vars);
    double user[n_user] = {};
    double error = 0;

    for (long i = 0; i < n_inner_cells; ++i) {
        sim->eqns->user.compute(x[i], sim->time, vars[i], user);
        error += volume[i] * (user[i_user] - vars[i][i_vars]) * (user[i_user] - vars[i][i_vars]);
    }
    error = sqrt(sync_sum(error) / total_volume);

    if (sim->eqns->mesh->rank == 0) {
        char key[128];
        snprintf(key, sizeof(key), "L2 error %s", sim->eqns->vars.name[i_vars]);
        printf(" | " FMT_KEY ": %g\n", key, error);
    }
}

static void terminate(int) { m_terminate = 1; }

static void write(Simulation *sim, double *output_time, long *output_iter)
{
    equations_write(sim->eqns, sim->prefix, sim->output_count, sim->time);
    sim->output_count += 1;
    *output_time = MAX(*output_time, sim->time + sim->output_time);
    *output_iter = MAX(*output_iter, sim->iter + sim->output_iter);
}
