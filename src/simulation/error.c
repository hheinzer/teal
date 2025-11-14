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

    scalar *norm = arena_malloc(stride, sizeof(*norm));
    equations_norm(sim->eqns, norm, time);
    scalar max_norm = array_fmax(norm, stride);

    int len = 0;
    for (number i = 0; i < num; i++) {
        len = lmax(len, strlen(name[i]));
    }

    println("Simulation error");
    for (number j = 0, i = 0; i < num; j += type[i++]) {
        switch (type[i]) {
            case SCALAR: println("\t %-*s : %g", len, name[i], norm[j + 0]); break;
            case VECTOR:
                println("\t %-*s : %g %g %g", len, name[i], norm[j + 0], norm[j + 1], norm[j + 2]);
                break;
            case MATRIX:
                println("\t %-*s : %g %g %g  %g %g %g  %g %g %g", len, name[i], norm[j + 0],
                        norm[j + 1], norm[j + 2], norm[j + 3], norm[j + 4], norm[j + 5],
                        norm[j + 6], norm[j + 7], norm[j + 8]);
                break;
            default: assert(false);
        }
    }

    arena_load(save);
    return max_norm;
}
