#include <assert.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <signal.h>

#include "simulation.h"
#include "teal/arena.h"
#include "teal/array.h"
#include "teal/sync.h"
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
    scalar courant = sim->advance.courant;

    scalar max_time = sim->time.max;
    scalar out_time = sim->time.out;

    long max_iter = sim->iter.max;
    long out_iter = sim->iter.out;

    const char *term_condition = sim->termination.condition;
    long variable = sim->termination.variable;
    scalar threshold = sim->termination.threshold;

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

    sync.wait = 0;
    double wtime_beg = MPI_Wtime();
    double wtime_last = wtime_beg;

    long iter = 0;
    while (iter < max_iter && is_less(time, max_time) && !has_converged && !sig_terminate) {
        scalar max_step = fmin(max_time, out_time) - time;
        scalar step0 = advance(eqns, &time, residual, max_step, courant, ctx);
        iter += 1;

        assert(isfinite(step0));
        for (long i = 0; i < len; i++) {
            assert(isfinite(residual[i]));
        }

        scalar max_residual = array_fmax(residual, len);
        if (term_condition) {
            if (variable >= 0) {
                max_residual = residual[variable];
            }
            has_converged = (max_residual < threshold);
        }

        bool has_reached_time = is_close_or_greater(time, fmin(max_time, out_time));
        bool has_reached_iter = iter >= lmin(max_iter, out_iter);
        if (has_reached_time || has_reached_iter || has_converged || sig_terminate) {
            double wtime_now = MPI_Wtime();
            double wtime = wtime_now - wtime_last;
            if (prefix) {
                equations_write(eqns, prefix, time, index++);
            }
            if (is_less(sim->time.out, SCALAR_MAX)) {
                out_time += sim->time.out;
            }
            if (sim->iter.out < LONG_MAX) {
                out_iter += sim->iter.out;
            }
            println("\t %13ld %13g %13g %13g %13g", iter, time, step0, max_residual, wtime);
            wtime_last = wtime_now;
        }
    }

    double wtime_end = MPI_Wtime();

    signal(SIGINT, SIG_DFL);
    signal(SIGTERM, SIG_DFL);

    char wtime[128];
    seconds_to_str(wtime, wtime_end - wtime_beg);
    println("\t computation time  : %s", wtime);

    char avg_wait[128];
    char max_wait[128];
    seconds_to_str(avg_wait, sync_fsum(sync.wait) / sync.size);
    seconds_to_str(max_wait, sync_fmax(sync.wait));
    println("\t avg/max wait time : %s / %s", avg_wait, max_wait);

    arena_load(save);
    return time;
}
