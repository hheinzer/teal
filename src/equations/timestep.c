#include <math.h>

#include "equations.h"
#include "teal/assert.h"
#include "teal/sync.h"

scalar equations_timestep(const Equations *eqns, const void *variable_, scalar *step)
{
    assert(eqns && variable_);

    number num = eqns->mesh->cells.num_inner;
    scalar *volume = eqns->mesh->cells.volume;
    vector *projection = eqns->mesh->cells.projection;

    number stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Timestep *timestep = eqns->timestep;

    const scalar(*variable)[stride] = variable_;

    scalar min = INFINITY;
    if (step) {
        for (number i = 0; i < num; i++) {
            step[i] = timestep(variable[i], property, volume[i], projection[i]);
            min = fmin(min, step[i]);
        }
    }
    else {
        for (number i = 0; i < num; i++) {
            min = fmin(min, timestep(variable[i], property, volume[i], projection[i]));
        }
    }
    return sync_fmin(min);
}
