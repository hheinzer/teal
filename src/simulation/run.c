#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <signal.h>

#include "equations.h"
#include "mesh.h"
#include "simulation.h"
#include "teal.h"
#include "utils.h"

static volatile sig_atomic_t sig_terminate = 0;

static void handler(int sig)
{
    (void)sig;
    sig_terminate = 1;
}

double simulation_run(Simulation *sim)
{
    assert(sim);

    const Equations *eqns = sim->eqns;
    int stride = eqns->conserved.stride;

    const char *prefix = sim->prefix;
    double courant = sim->advance.courant;
    void *context = sim->advance.context;
    Advance *advance = sim->advance.method;

    double max_time = sim->time.max;
    double out_time = sim->time.out;
    int max_iter = sim->iter.max;
    int out_iter = sim->iter.out;

    const char *condition = sim->termination.condition;
    int variable = sim->termination.variable;
    double threshold = sim->termination.threshold;

    double time = 0;
    int index = 0;

    if (prefix[0]) {
        mesh_write(eqns->mesh, prefix);
        equations_write(eqns, prefix, time, index++);
    }
    out_time += time;

    double *residual = teal_calloc(stride, sizeof(*residual));

    teal_print("Running simulation");
    teal_print("\t %13s %13s %13s %13s %13s", "iter", "time", "step", "residual", "wtime");

    signal(SIGINT, handler);
    signal(SIGTERM, handler);

    double wtime_beg = MPI_Wtime();
    double wtime_last = wtime_beg;

    int iter = 0;
    int has_converged = 0;
    while (iter < max_iter && is_less(time, max_time) && !has_converged && !sig_terminate) {
        double max_step = fmin(max_time, out_time) - time;
        double step0 = advance(eqns, &time, residual, max_step, courant, context);
        iter += 1;

        double max_residual = 0;
        for (int i = 0; i < stride; i++) {
            if (residual[i] > max_residual) {
                max_residual = residual[i];
            }
        }

        if (condition[0]) {
            has_converged = ((variable >= 0) ? residual[variable] : max_residual) < threshold;
        }

        int has_reached_time = is_greater_equal(time, fmin(max_time, out_time));
        int has_reached_iter = iter >= (out_iter < max_iter ? out_iter : max_iter);

        if (has_reached_time || has_reached_iter || has_converged || sig_terminate) {
            double wtime_now = MPI_Wtime();
            if (prefix[0]) {
                equations_write(eqns, prefix, time, index++);
            }
            if (sim->time.out < DBL_MAX) {
                out_time += sim->time.out;
            }
            if (sim->iter.out < INT_MAX) {
                out_iter += sim->iter.out;
            }
            teal_print("\t %13d %13g %13g %13g %13g", iter, time, step0, max_residual,
                       wtime_now - wtime_last);
            wtime_last = wtime_now;
        }
    }

    double wtime_end = MPI_Wtime();

    signal(SIGINT, SIG_DFL);
    signal(SIGTERM, SIG_DFL);

    teal_print("\t computation time : %g", wtime_end - wtime_beg);

    teal_free(residual);
    return time;
}
