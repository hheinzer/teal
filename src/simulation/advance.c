#include "advance.h"

#include <math.h>

#include "teal/arena.h"
#include "teal/utils.h"

scalar explicit_euler(const Equations *eqns, scalar *time, void *residual_, scalar courant,
                      scalar max_timestep, void *context)
{
    unused(context);

    Arena save = arena_save();

    number num = eqns->mesh->cells.num_inner;
    number len = eqns->variables.len;
    number stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Update *primitive = eqns->variables.primitive;

    scalar(*variable)[stride] = eqns->variables.data;
    scalar(*derivative)[len] = equations_derivative(eqns, variable, *time);

    scalar full_timestep = courant * equations_timestep(eqns, variable, 0);
    scalar timestep = fmin(full_timestep, max_timestep);

    for (number i = 0; i < num; i++) {
        for (number j = 0; j < len; j++) {
            variable[i][j] += derivative[i][j] * timestep;
        }
        primitive(variable[i], property);
    }

    *time += timestep;
    equations_residual(eqns, derivative, residual_);

    arena_load(save);
    return full_timestep;
}
