#include <string.h>

#include "simulation.h"
#include "teal/arena.h"
#include "teal/array.h"
#include "teal/assert.h"
#include "teal/utils.h"

scalar simulation_error(const Simulation *sim, scalar time)
{
    assert(sim);

    Arena save = arena_save();

    number num = sim->eqns->variables.num;
    number stride = sim->eqns->variables.stride;
    Name *name = sim->eqns->variables.name;
    Type *type = sim->eqns->variables.type;

    scalar *error = arena_malloc(stride, sizeof(*error));
    equations_error(sim->eqns, error, time);
    scalar max_error = array_fmax(error, stride);

    int len = 0;
    for (number i = 0; i < num; i++) {
        len = lmax(len, strlen(name[i]));
    }

    println("Simulation error");
    for (number j = 0, i = 0; i < num; j += type[i++]) {
        switch (type[i]) {
            case SCALAR: println("\t %-*s : %g", len, name[i], error[j + 0]); break;
            case VECTOR:
                println("\t %-*s : %g %g %g", len, name[i], error[j + 0], error[j + 1],
                        error[j + 2]);
                break;
            case MATRIX:
                println("\t %-*s : %g %g %g  %g %g %g  %g %g %g", len, name[i], error[j + 0],
                        error[j + 1], error[j + 2], error[j + 3], error[j + 4], error[j + 5],
                        error[j + 6], error[j + 7], error[j + 8]);
                break;
            default: assert(false);
        }
    }

    arena_load(save);
    return max_error;
}
