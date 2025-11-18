#include <math.h>
#include <string.h>

#include "simulation.h"
#include "teal/assert.h"
#include "teal/utils.h"

void simulation_summary(const Simulation *sim)
{
    assert(sim);

    println("Simulation summary");
    if (sim->prefix) {
        println("\t prefix           : %s", sim->prefix);
    }
    println("\t courant          : %g", sim->courant);

    if (isfinite(sim->time.max)) {
        println("\t max time         : %g", sim->time.max);
    }
    if (isfinite(sim->time.output)) {
        println("\t output time      : %g", sim->time.output);
    }

    if (sim->iter.max < NUMBER_MAX) {
        println("\t max iter         : %td", sim->iter.max);
    }
    if (sim->iter.output < NUMBER_MAX) {
        println("\t output iter      : %td", sim->iter.output);
    }

    if (sim->termination.condition) {
        println("\t termination      : %s < %g", sim->termination.condition,
                sim->termination.residual);
    }

    println("\t advance          : %s", sim->advance.name);
    if (!strcmp(sim->advance.name, "lserk")) {
        const RungeKutta *ctx = sim->advance.ctx;
        println("\t time order       : %td", ctx->time_order);
        println("\t number of stages : %td", ctx->num_stages);
    }
    if (!strcmp(sim->advance.name, "implicit euler")) {
        const NewtonKrylov *ctx = sim->advance.ctx;
        println("\t newton tolerance : %g", ctx->newton_tolerance);
        println("\t krylov tolerance : %g", ctx->krylov_tolerance);
        println("\t krylov dimension : %td", ctx->krylov_dimension);
    }

    scalar step0 = sim->courant * equations_timestep(sim->eqns, sim->eqns->variables.data, 0);
    println("\t initial timestep : %g", step0);
}
