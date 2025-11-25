#include <assert.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <signal.h>

#include "simulation.h"
#include "teal/arena.h"
#include "teal/array.h"
#include "teal/utils.h"

static volatile sig_atomic_t sig_terminate = 0;

static void handler(int sig)
{
    (void)sig;
    sig_terminate = 1;
}

scalar simulation_run(Simulation *sim)
{
    assert(sim);
    Arena save = arena_save();

    const Equations *eqns = sim->eqns;
    long len = eqns->variables.len;

    const char *prefix = sim->prefix;
    scalar courant = sim->courant;

    scalar max_time = sim->time.max;
    scalar out_time = sim->time.output;

    long max_iter = sim->iter.max;
    long out_iter = sim->iter.output;

    const char *term_condition = sim->termination.condition;
    long term_variable = sim->termination.variable;
    scalar term_residual = sim->termination.residual;

    const void *ctx = sim->advance.ctx;
    Advance *advance = sim->advance.method;

    scalar time;
    long index;
    equations_restart(eqns, &time, &index);
    if (prefix) {
        mesh_write(eqns->mesh, prefix);
        equations_write(eqns, prefix, time, index++);
    }
    out_time += time;

    bool has_converged = false;
    scalar *residual = arena_malloc(len, sizeof(*residual));

    println("Running simulation");
    println("\t %13s %13s %13s %13s %13s", "iter", "time", "step", "residual", "wtime");

    signal(SIGINT, handler);
    signal(SIGTERM, handler);

    scalar wtime_beg = MPI_Wtime();
    scalar wtime_last = wtime_beg;

    for (long iter = 0; iter < max_iter && time < max_time && !has_converged && !sig_terminate;) {
        scalar max_step = fmin(max_time, out_time) - time;
        scalar step0 = advance(eqns, &time, residual, courant, max_step, ctx);

        if (!isfinite(step0)) {
            error("non-finite timestep (%ld)", iter);
        }

        for (long i = 0; i < len; i++) {
            if (!isfinite(residual[i])) {
                error("non-finite residual (%ld)", iter);
            }
        }

        iter += 1;

        scalar max_residual = array_fmax(residual, len);
        if (term_condition) {
            if (term_variable >= 0) {
                max_residual = residual[term_variable];
            }
            has_converged = (max_residual < term_residual);
        }

        if (time >= fmin(max_time, out_time) || iter >= lmin(max_iter, out_iter) || has_converged ||
            sig_terminate) {
            scalar wtime_now = MPI_Wtime();
            scalar wtime = wtime_now - wtime_last;
            if (prefix) {
                equations_write(eqns, prefix, time, index++);
            }
            if (isfinite(sim->time.output)) {
                out_time = time + sim->time.output;
            }
            if (sim->iter.output < LONG_MAX) {
                out_iter = iter + sim->iter.output;
            }
            println("\t %13ld %13g %13g %13g %13g", iter, time, step0, max_residual, wtime);
            wtime_last = wtime_now;
        }
    }

    scalar wtime_end = MPI_Wtime();

    signal(SIGINT, SIG_DFL);
    signal(SIGTERM, SIG_DFL);

    char wtime[128];
    seconds_to_str(wtime, wtime_end - wtime_beg);
    println("\t computation time : %s", wtime);

    arena_load(save);
    return time;
}
