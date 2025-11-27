#include <assert.h>
#include <limits.h>
#include <math.h>
#include <string.h>

#include "simulation.h"
#include "teal/utils.h"

void simulation_summary(const Simulation *sim)
{
    assert(sim);

    println("Simulation summary");
    println("\t prefix            : %s", sim->prefix);
    println("\t courant           : %g", sim->courant);

    if (isfinite(sim->time.max)) {
        println("\t max time          : %g", sim->time.max);
    }
    if (isfinite(sim->time.out)) {
        println("\t output time       : %g", sim->time.out);
    }

    if (sim->iter.max < LONG_MAX) {
        println("\t max iter          : %ld", sim->iter.max);
    }
    if (sim->iter.out < LONG_MAX) {
        println("\t output iter       : %ld", sim->iter.out);
    }

    if (sim->termination.condition) {
        println("\t termination       : %s < %g", sim->termination.condition,
                sim->termination.residual);
    }

    println("\t advance method    : %s", sim->advance.name);
    if (!strcmp(sim->advance.name, "lserk")) {
        const RungeKutta *ctx = sim->advance.ctx;
        println("\t time order        : %ld", ctx->time_order);
        println("\t number of stages  : %ld", ctx->num_stages);
    }
    if (!strcmp(sim->advance.name, "implicit euler")) {
        const NewtonKrylov *ctx = sim->advance.ctx;
        println("\t newton tolerance  : %g", ctx->newton_tolerance);
        println("\t krylov tolerance  : %g", ctx->krylov_tolerance);
        println("\t krylov dimension  : %ld", ctx->krylov_dimension);
    }

    scalar step0 = sim->courant * equations_time_step(sim->eqns, sim->eqns->variables.data, 0);
    println("\t initial time step : %g", step0);
}
