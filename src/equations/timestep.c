#include <assert.h>
#include <math.h>

#include "equations.h"
#include "teal/sync.h"

scalar equations_timestep(const Equations *eqns, const void *variable_, scalar *step)
{
    assert(eqns && variable_);

    long num_inner = eqns->mesh->cells.num_inner;
    scalar *volume = eqns->mesh->cells.volume;
    vector *projection = eqns->mesh->cells.projection;

    long stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Timestep *timestep = eqns->timestep;

    const scalar(*variable)[stride] = variable_;

    scalar min = INFINITY;
    if (step) {
        for (long i = 0; i < num_inner; i++) {
            step[i] = timestep(variable[i], property, volume[i], projection[i]);
            min = fmin(min, step[i]);
        }
    }
    else {
        for (long i = 0; i < num_inner; i++) {
            min = fmin(min, timestep(variable[i], property, volume[i], projection[i]));
        }
    }
    return sync_fmin(min);
}
