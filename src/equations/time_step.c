#include <assert.h>
#include <math.h>

#include "equations.h"
#include "teal/sync.h"
#include "teal/utils.h"

scalar equations_time_step(const Equations *eqns, const void *variable_)
{
    assert(eqns && variable_);

    long num_inner = eqns->mesh->cells.num_inner;
    scalar *volume = eqns->mesh->cells.volume;
    vector *projection = eqns->mesh->cells.projection;

    long stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    TimeStep *time_step = eqns->time_step;

    const scalar(*variable)[stride] = variable_;

    scalar min = SCALAR_MAX;
    for (long i = 0; i < num_inner; i++) {
        min = fmin(min, time_step(variable[i], property, volume[i], projection[i]));
    }
    return sync_fmin(min);
}
