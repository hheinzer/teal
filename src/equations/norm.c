#include <assert.h>
#include <math.h>

#include "equations.h"
#include "teal/arena.h"
#include "teal/sync.h"
#include "teal/utils.h"

void *equations_norm(const Equations *eqns, scalar time)
{
    assert(eqns);

    Arena save = arena_save();

    long num = eqns->mesh->cells.num_inner;
    vector *center = eqns->mesh->cells.center;
    scalar *volume = eqns->mesh->cells.volume;
    scalar sum_volume = eqns->mesh->cells.sum_volume;

    long stride = eqns->variables.stride;
    scalar(*variable)[stride] = eqns->variables.data;
    scalar *property = eqns->properties.data;
    Compute *compute = eqns->user.compute;
    Update *conserved = eqns->user.conserved;

    scalar *norm = arena_calloc(stride, sizeof(*norm));
    scalar *user = arena_calloc(stride, sizeof(*user));

    for (long i = 0; i < num; i++) {
        compute(user, property, center[i], time, variable[i]);
        conserved(user, property);
        for (long j = 0; j < stride; j++) {
            norm[j] += volume[i] * sq(user[j] - variable[i][j]);
        }
    }
    for (long i = 0; i < stride; i++) {
        norm[i] = sqrt(sync_fsum(norm[i]) / sum_volume);
    }

    arena_load(save);
    return arena_smuggle(norm, stride, sizeof(*norm));
}
