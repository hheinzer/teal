#include <math.h>
#include <stdint.h>
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

    if (sim->iter.max < PTRDIFF_MAX) {
        println("\t max iter         : %td", sim->iter.max);
    }
    if (sim->iter.output < PTRDIFF_MAX) {
        println("\t output iter      : %td", sim->iter.output);
    }

    if (sim->termination.condition) {
        println("\t termination      : %s < %g", sim->termination.condition,
                sim->termination.residual);
    }

    println("\t advance          : %s", sim->advance.name);

    scalar timestep = equations_timestep(sim->eqns, sim->eqns->variables.data, 0);
    println("\t initial timestep : %g", sim->courant * timestep);
}
