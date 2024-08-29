#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>

#include "simulation.h"
#include "teal/array.h"
#include "teal/option.h"
#include "teal/print.h"
#include "teal/sync.h"
#include "teal/utils.h"

static volatile sig_atomic_t m_terminate = 0;

static void terminate(int);

static void print_header(void);

static void print(const Simulation *sim, double dt, double residual, double wtime);

static void write(Simulation *sim, double *output_time, long *output_iter);

void simulation_run(Simulation *sim)
{
    double output_time = sim->output_time;
    long output_iter = sim->output_iter;
    mesh_write(sim->eqns->mesh, sim->prefix);
    write(sim, &output_time, &output_iter);

    const long n_cons = sim->eqns->n_cons;
    const long abort = sim->abort_variable;
    const double abort_residual = sim->abort_residual;
    double residual[n_cons];
    int converged = 0;

    if (sync.rank == 0) print_header();

    signal(SIGINT, terminate);
    signal(SIGTERM, terminate);

    const double timer_start = MPI_Wtime();
    double timer_last = timer_start;

    while (sim->time < sim->max_time && sim->iter < sim->max_iter && !converged && !m_terminate) {
        const double max_dt = min(sim->max_time, output_time) - sim->time;
        const double dt = sim->advance.method(sim, max_dt);
        assert(isfinite(dt));
        sim->iter += 1;

        equations_residual(sim->eqns, residual);
        for (long i = 0; i < n_cons; ++i) assert(isfinite(residual[i]));
        const double max_residual = (abort < 0 ? array_max(residual, n_cons) : residual[abort]);
        converged = (max_residual < abort_residual);

        if (sim->time >= min(sim->max_time, output_time) ||
            sim->iter >= min(sim->max_iter, output_iter) || converged || m_terminate) {
            const double timer_now = MPI_Wtime();
            write(sim, &output_time, &output_iter);
            if (sync.rank == 0) print(sim, dt, max_residual, timer_now - timer_last);
            timer_last = timer_now;
        }
    }

    const double timer_stop = MPI_Wtime();

    signal(SIGINT, SIG_DFL);
    signal(SIGTERM, SIG_DFL);

    if (sync.rank == 0 && !option.quiet) {
        print_key("computation time", "%g", timer_stop - timer_start);
        if (sim->iter_newton) print_key("Newton iterations", "%ld", sim->iter_newton);
        if (sim->iter_krylov) print_key("Krylov iterations", "%ld", sim->iter_krylov);
    }
}

static void terminate(int)
{
    m_terminate = 1;
}

static void print_header(void)
{
    if (option.quiet) return;
    printf(" |  %13s  %13s  %13s  %13s  %13s\n", "iter", "time", "dt", "residual", "wtime");
}

static void print(const Simulation *sim, double dt, double residual, double wtime)
{
    if (option.quiet) return;
    printf(" |  %13ld  %13g  %13g  %13g  %13g\n", sim->iter, sim->time, dt, residual, wtime);
}

static void write(Simulation *sim, double *output_time, long *output_iter)
{
    equations_write(sim->eqns, sim->prefix, sim->output_count, sim->time);
    sim->output_count += 1;
    if (sim->output_time < DBL_MAX) *output_time = sim->time + sim->output_time;
    if (sim->output_iter < LONG_MAX) *output_iter = sim->iter + sim->output_iter;
}
