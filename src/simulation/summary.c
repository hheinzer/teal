#include <assert.h>
#include <float.h>
#include <limits.h>

#include "advance.h"
#include "simulation.h"
#include "teal.h"

void simulation_summary(const Simulation *sim)
{
    assert(sim);

    teal_print("Simulation summary");
    teal_print("\t prefix           : %s", sim->prefix);

    if (sim->time.max < DBL_MAX) {
        teal_print("\t max time         : %g", sim->time.max);
    }
    if (sim->time.out < DBL_MAX) {
        teal_print("\t out time         : %g", sim->time.out);
    }

    if (sim->iter.max < INT_MAX) {
        teal_print("\t max iter         : %d", sim->iter.max);
    }
    if (sim->iter.out < INT_MAX) {
        teal_print("\t out iter         : %d", sim->iter.out);
    }

    if (sim->termination.condition[0]) {
        teal_print("\t termination      : %s < %g", sim->termination.condition,
                   sim->termination.threshold);
    }

    teal_print("\t advance method   : %s", sim->advance.name);
    teal_print("\t courant number   : %g", sim->advance.courant);
    if (sim->advance.method == lserk) {
        RungeKutta *context = sim->advance.context;
        teal_print("\t time order       : %d", context->time_order);
        teal_print("\t number of stages : %d", context->num_stages);
    }
    if (sim->advance.method == implicit_euler) {
        NewtonKrylov *context = sim->advance.context;
        teal_print("\t newton tolerance : %g", context->newton_tolerance);
        teal_print("\t krylov tolerance : %g", context->krylov_tolerance);
        teal_print("\t krylov dimension : %d", context->krylov_dimension);
    }
}
