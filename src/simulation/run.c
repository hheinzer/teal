#include <math.h>
#include <signal.h>
#include <stdio.h>

#include "core/array.h"
#include "core/utils.h"
#include "simulation.h"
#include "teal.h"

static volatile sig_atomic_t m_terminate = 0;

static void terminate(int signum);

static void print_header(void);

static void print(const Simulation *sim, double dt, double residual, double wtime);

static void write(Simulation *sim, double *output_time, long *output_iter);

void simulation_run(Simulation *sim)
{
    double output_time = sim->output_time;
    long output_iter = sim->output_iter;
    mesh_write(sim->eqns->mesh, sim->prefix);
    write(sim, &output_time, &output_iter);

    const long n_vars = sim->eqns->n_vars;
    const long a_var = sim->abort_variable;
    const double a_res = sim->abort_residual;
    double residual[n_vars];
    int converged = 0;

    if (teal.rank == 0) print_header();

    signal(SIGINT, terminate);
    signal(SIGTERM, terminate);

    const double timer_start = MPI_Wtime();
    double timer_last = timer_start;

    while (sim->time < sim->max_time && sim->iter < sim->max_iter && !converged && !m_terminate) {
        const double max_dt = min(sim->max_time, output_time) - sim->time;
        const double dt = sim->advance.func(sim, max_dt);
        ensure(isfinite(dt));
        sim->iter += 1;

        equations_residual(sim->eqns, residual);
        const double max_residual = (a_var < 0 ? array_max(residual, n_vars) : residual[a_var]);
        converged = (max_residual < a_res);

        if (sim->time >= min(sim->max_time, output_time) ||
            sim->iter >= min(sim->max_iter, output_iter) || converged || m_terminate) {
            const double timer_now = MPI_Wtime();
            write(sim, &output_time, &output_iter);
            if (teal.rank == 0) print(sim, dt, max_residual, timer_now - timer_last);
            timer_last = timer_now;
        }
    }

    const double timer_stop = MPI_Wtime();

    signal(SIGINT, SIG_DFL);
    signal(SIGTERM, SIG_DFL);

    if (teal.rank == 0 && !teal.quiet)
        printf(" | " KEYFMT ": %g\n", "computation time", timer_stop - timer_start);
}

static void terminate(int)
{
    m_terminate = 1;
}

static void print_header(void)
{
    if (teal.quiet) return;
    printf(" |  %13s  %13s  %13s  %13s  %13s\n", "iter", "time", "dt", "residual", "wtime");
}

static void print(const Simulation *sim, double dt, double residual, double wtime)
{
    if (teal.quiet) return;
    printf(" |  %13ld  %13g  %13g  %13g  %13g\n", sim->iter, sim->time, dt, residual, wtime);
}

static void write(Simulation *sim, double *output_time, long *output_iter)
{
    equations_write(sim->eqns, sim->prefix, sim->output_count, sim->time, sim->iter);
    sim->output_count += 1;
    *output_time = max(*output_time, sim->time + sim->output_time);
    *output_iter = max(*output_iter, sim->iter + sim->output_iter);
}
