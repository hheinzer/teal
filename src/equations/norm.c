#include <math.h>
#include <string.h>

#include "equations.h"
#include "teal/arena.h"
#include "teal/assert.h"
#include "teal/sync.h"
#include "teal/utils.h"

void equations_norm(const Equations *eqns, void *norm_, scalar time)
{
    assert(eqns && norm_);

    Arena save = arena_save();

    number num = eqns->mesh->cells.num_inner;
    vector *center = eqns->mesh->cells.center;
    scalar *volume = eqns->mesh->cells.volume;
    scalar sum_volume = eqns->mesh->cells.sum_volume;

    number stride = eqns->variables.stride;
    scalar(*variable)[stride] = eqns->variables.data;
    scalar *property = eqns->properties.data;
    Compute *compute = eqns->user.compute;
    Update *conserved = eqns->user.conserved;

    scalar *user = arena_calloc(stride, sizeof(*user));
    scalar *norm = norm_;

    memset(norm, 0, stride * sizeof(*norm));
    for (number i = 0; i < num; i++) {
        compute(user, variable[i], property, center[i], time);
        conserved(user, property);
        for (number j = 0; j < stride; j++) {
            norm[j] += volume[i] * pow2(user[j] - variable[i][j]);
        }
    }
    for (number i = 0; i < stride; i++) {
        norm[i] = sqrt(sync_fsum(norm[i]) / sum_volume);
    }

    arena_load(save);
}
