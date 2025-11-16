#include "advance.h"

#include <math.h>

#include "teal/arena.h"
#include "teal/utils.h"

scalar euler(const Equations *eqns, scalar *time, void *residual_, scalar courant,
             scalar max_timestep, const void *context_)
{
    unused(context_);
    Arena save = arena_save();

    number num = eqns->mesh->cells.num_inner;
    number len = eqns->variables.len;
    number stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Update *primitive = eqns->variables.primitive;

    scalar(*variable)[stride] = eqns->variables.data;
    scalar(*derivative)[len] = equations_derivative(eqns, variable, 0, *time);

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

scalar midpoint(const Equations *eqns, scalar *time, void *residual_, scalar courant,
                scalar max_timestep, const void *context_)
{
    unused(context_);
    Arena save = arena_save();

    number num_cells = eqns->mesh->cells.num;
    number num = eqns->mesh->cells.num_inner;
    number len = eqns->variables.len;
    number stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Update *primitive = eqns->variables.primitive;

    scalar(*variable)[stride] = eqns->variables.data;
    scalar(*derivative)[len] = equations_derivative(eqns, variable, 0, *time);

    scalar full_timestep = courant * equations_timestep(eqns, variable, 0);
    scalar timestep = fmin(full_timestep, max_timestep);

    scalar(*variable1)[stride] = arena_malloc(num_cells, sizeof(*variable1));
    for (number i = 0; i < num; i++) {
        for (number j = 0; j < len; j++) {
            variable1[i][j] = variable[i][j] + (derivative[i][j] * timestep / 2);
        }
        primitive(variable1[i], property);
    }
    scalar(*derivative1)[len] = equations_derivative(eqns, variable1, 0, *time + (timestep / 2));

    for (number i = 0; i < num; i++) {
        for (number j = 0; j < len; j++) {
            variable[i][j] += derivative1[i][j] * timestep;
        }
        primitive(variable[i], property);
    }

    *time += timestep;
    equations_residual(eqns, derivative, residual_);

    arena_load(save);
    return full_timestep;
}

scalar heun(const Equations *eqns, scalar *time, void *residual_, scalar courant,
            scalar max_timestep, const void *context_)
{
    unused(context_);
    Arena save = arena_save();

    number num_cells = eqns->mesh->cells.num;
    number num = eqns->mesh->cells.num_inner;
    number len = eqns->variables.len;
    number stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Update *primitive = eqns->variables.primitive;

    scalar(*variable)[stride] = eqns->variables.data;
    scalar(*derivative)[len] = equations_derivative(eqns, variable, 0, *time);

    scalar full_timestep = courant * equations_timestep(eqns, variable, 0);
    scalar timestep = fmin(full_timestep, max_timestep);

    scalar(*variable1)[stride] = arena_malloc(num_cells, sizeof(*variable1));
    for (number i = 0; i < num; i++) {
        for (number j = 0; j < len; j++) {
            variable1[i][j] = variable[i][j] + (derivative[i][j] * timestep);
        }
        primitive(variable1[i], property);
    }
    scalar(*derivative1)[len] = equations_derivative(eqns, variable1, 0, *time + timestep);

    for (number i = 0; i < num; i++) {
        for (number j = 0; j < len; j++) {
            variable[i][j] += (derivative[i][j] + derivative1[i][j]) * timestep / 2;
        }
        primitive(variable[i], property);
    }

    *time += timestep;
    equations_residual(eqns, derivative, residual_);

    arena_load(save);
    return full_timestep;
}

scalar ralston(const Equations *eqns, scalar *time, void *residual_, scalar courant,
               scalar max_timestep, const void *context_)
{
    unused(context_);
    Arena save = arena_save();

    number num_cells = eqns->mesh->cells.num;
    number num = eqns->mesh->cells.num_inner;
    number len = eqns->variables.len;
    number stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Update *primitive = eqns->variables.primitive;

    scalar(*variable)[stride] = eqns->variables.data;
    scalar(*derivative)[len] = equations_derivative(eqns, variable, 0, *time);

    scalar full_timestep = courant * equations_timestep(eqns, variable, 0);
    scalar timestep = fmin(full_timestep, max_timestep);

    scalar(*variable1)[stride] = arena_malloc(num_cells, sizeof(*variable1));
    for (number i = 0; i < num; i++) {
        for (number j = 0; j < len; j++) {
            variable1[i][j] = variable[i][j] + (derivative[i][j] * 2 * timestep / 3);
        }
        primitive(variable1[i], property);
    }
    scalar(*derivative1)[len] =
        equations_derivative(eqns, variable1, 0, *time + (2 * timestep / 3));

    for (number i = 0; i < num; i++) {
        for (number j = 0; j < len; j++) {
            variable[i][j] += timestep * (derivative[i][j] + 3 * derivative1[i][j]) / 4;
        }
        primitive(variable[i], property);
    }

    *time += timestep;
    equations_residual(eqns, derivative, residual_);

    arena_load(save);
    return full_timestep;
}

scalar ssprk2(const Equations *eqns, scalar *time, void *residual_, scalar courant,
              scalar max_timestep, const void *context_)
{
    unused(context_);
    Arena save = arena_save();

    number num_cells = eqns->mesh->cells.num;
    number num = eqns->mesh->cells.num_inner;
    number len = eqns->variables.len;
    number stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Update *primitive = eqns->variables.primitive;

    scalar(*variable)[stride] = eqns->variables.data;
    scalar(*derivative)[len] = equations_derivative(eqns, variable, 0, *time);

    scalar full_timestep = courant * equations_timestep(eqns, variable, 0);
    scalar timestep = fmin(full_timestep, max_timestep);

    scalar(*variable1)[stride] = arena_malloc(num_cells, sizeof(*variable1));
    for (number i = 0; i < num; i++) {
        for (number j = 0; j < len; j++) {
            variable1[i][j] = variable[i][j] + (derivative[i][j] * timestep);
        }
        primitive(variable1[i], property);
    }
    scalar(*derivative1)[len] = equations_derivative(eqns, variable1, 0, *time + timestep);

    for (number i = 0; i < num; i++) {
        for (number j = 0; j < len; j++) {
            variable[i][j] = (variable[i][j] + variable1[i][j] + timestep * derivative1[i][j]) / 2;
        }
        primitive(variable[i], property);
    }

    *time += timestep;
    equations_residual(eqns, derivative, residual_);

    arena_load(save);
    return full_timestep;
}

scalar ssprk3(const Equations *eqns, scalar *time, void *residual_, scalar courant,
              scalar max_timestep, const void *context_)
{
    unused(context_);
    Arena save = arena_save();

    number num_cells = eqns->mesh->cells.num;
    number num = eqns->mesh->cells.num_inner;
    number len = eqns->variables.len;
    number stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Update *primitive = eqns->variables.primitive;

    scalar(*variable)[stride] = eqns->variables.data;
    scalar(*derivative)[len] = equations_derivative(eqns, variable, 0, *time);

    scalar full_timestep = courant * equations_timestep(eqns, variable, 0);
    scalar timestep = fmin(full_timestep, max_timestep);

    scalar(*variable1)[stride] = arena_malloc(num_cells, sizeof(*variable1));
    for (number i = 0; i < num; i++) {
        for (number j = 0; j < len; j++) {
            variable1[i][j] = variable[i][j] + (derivative[i][j] * timestep);
        }
        primitive(variable1[i], property);
    }
    scalar(*derivative1)[len] = equations_derivative(eqns, variable1, 0, *time + timestep);

    scalar(*variable2)[stride] = arena_malloc(num_cells, sizeof(*variable2));
    for (number i = 0; i < num; i++) {
        for (number j = 0; j < len; j++) {
            variable2[i][j] =
                (3 * variable[i][j] + variable1[i][j] + timestep * derivative1[i][j]) / 4;
        }
        primitive(variable2[i], property);
    }
    scalar(*derivative2)[len] = equations_derivative(eqns, variable2, 0, *time + (timestep / 2));

    for (number i = 0; i < num; i++) {
        for (number j = 0; j < len; j++) {
            variable[i][j] =
                (variable[i][j] + 2 * (variable2[i][j] + timestep * derivative2[i][j])) / 3;
        }
        primitive(variable[i], property);
    }

    *time += timestep;
    equations_residual(eqns, derivative, residual_);

    arena_load(save);
    return full_timestep;
}

scalar rk3(const Equations *eqns, scalar *time, void *residual_, scalar courant,
           scalar max_timestep, const void *context_)
{
    unused(context_);
    Arena save = arena_save();

    number num_cells = eqns->mesh->cells.num;
    number num = eqns->mesh->cells.num_inner;
    number len = eqns->variables.len;
    number stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Update *primitive = eqns->variables.primitive;

    scalar(*variable)[stride] = eqns->variables.data;
    scalar(*derivative)[len] = equations_derivative(eqns, variable, 0, *time);

    scalar full_timestep = courant * equations_timestep(eqns, variable, 0);
    scalar timestep = fmin(full_timestep, max_timestep);

    scalar(*variable1)[stride] = arena_malloc(num_cells, sizeof(*variable1));
    for (number i = 0; i < num; i++) {
        for (number j = 0; j < len; j++) {
            variable1[i][j] = variable[i][j] + (derivative[i][j] * timestep / 2);
        }
        primitive(variable1[i], property);
    }
    scalar(*derivative1)[len] = equations_derivative(eqns, variable1, 0, *time + (timestep / 2));

    scalar(*variable2)[stride] = arena_malloc(num_cells, sizeof(*variable2));
    for (number i = 0; i < num; i++) {
        for (number j = 0; j < len; j++) {
            variable2[i][j] =
                variable[i][j] - (derivative[i][j] * timestep) + (derivative1[i][j] * 2 * timestep);
        }
        primitive(variable2[i], property);
    }
    scalar(*derivative2)[len] = equations_derivative(eqns, variable2, 0, *time + timestep);

    for (number i = 0; i < num; i++) {
        for (number j = 0; j < len; j++) {
            variable[i][j] +=
                timestep * (derivative[i][j] + 4 * derivative1[i][j] + derivative2[i][j]) / 6;
        }
        primitive(variable[i], property);
    }

    *time += timestep;
    equations_residual(eqns, derivative, residual_);

    arena_load(save);
    return full_timestep;
}

scalar rk4(const Equations *eqns, scalar *time, void *residual_, scalar courant,
           scalar max_timestep, const void *context_)
{
    unused(context_);
    Arena save = arena_save();

    number num_cells = eqns->mesh->cells.num;
    number num = eqns->mesh->cells.num_inner;
    number len = eqns->variables.len;
    number stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Update *primitive = eqns->variables.primitive;

    scalar(*variable)[stride] = eqns->variables.data;
    scalar(*derivative)[len] = equations_derivative(eqns, variable, 0, *time);

    scalar full_timestep = courant * equations_timestep(eqns, variable, 0);
    scalar timestep = fmin(full_timestep, max_timestep);

    scalar(*variable1)[stride] = arena_malloc(num_cells, sizeof(*variable1));
    for (number i = 0; i < num; i++) {
        for (number j = 0; j < len; j++) {
            variable1[i][j] = variable[i][j] + (derivative[i][j] * timestep / 2);
        }
        primitive(variable1[i], property);
    }
    scalar(*derivative1)[len] = equations_derivative(eqns, variable1, 0, *time + (timestep / 2));

    scalar(*variable2)[stride] = arena_malloc(num_cells, sizeof(*variable2));
    for (number i = 0; i < num; i++) {
        for (number j = 0; j < len; j++) {
            variable2[i][j] = variable[i][j] + (derivative1[i][j] * timestep / 2);
        }
        primitive(variable2[i], property);
    }
    scalar(*derivative2)[len] = equations_derivative(eqns, variable2, 0, *time + (timestep / 2));

    scalar(*variable3)[stride] = arena_malloc(num_cells, sizeof(*variable3));
    for (number i = 0; i < num; i++) {
        for (number j = 0; j < len; j++) {
            variable3[i][j] = variable[i][j] + (derivative2[i][j] * timestep);
        }
        primitive(variable3[i], property);
    }
    scalar(*derivative3)[len] = equations_derivative(eqns, variable3, 0, *time + timestep);

    for (number i = 0; i < num; i++) {
        for (number j = 0; j < len; j++) {
            variable[i][j] += timestep *
                              (derivative[i][j] + 2 * (derivative1[i][j] + derivative2[i][j]) +
                               derivative3[i][j]) /
                              6;
        }
        primitive(variable[i], property);
    }

    *time += timestep;
    equations_residual(eqns, derivative, residual_);

    arena_load(save);
    return full_timestep;
}

scalar lserk(const Equations *eqns, scalar *time, void *residual_, scalar courant,
             scalar max_timestep, const void *context_)
{
    unused(context_);
    Arena save = arena_save();

    static const scalar alpha_[3][6][7] = {
        {
            {0.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.3333, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.1481, 0.4000, 1.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.0833, 0.2069, 0.4265, 1.0000, 0.0000, 0.0000},
            {0.0000, 0.0533, 0.1263, 0.2375, 0.4414, 1.0000, 0.0000},
            {0.0000, 0.0370, 0.0851, 0.1521, 0.2562, 0.4512, 1.0000},
        },
        {
            {0.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},  // first order
            {0.0000, 0.4242, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.1918, 0.4929, 1.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.1084, 0.2602, 0.5052, 1.0000, 0.0000, 0.0000},
            {0.0000, 0.0695, 0.1602, 0.2898, 0.5060, 1.0000, 0.0000},
            {0.0000, 0.0482, 0.1085, 0.1885, 0.3050, 0.5063, 1.0000},
        },
        {
            {0.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},  // first order
            {0.0000, 0.6612, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.2884, 0.5010, 1.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.1666, 0.3027, 0.5275, 1.0000, 0.0000, 0.0000},
            {0.0000, 0.1067, 0.1979, 0.3232, 0.5201, 1.0000, 0.0000},
            {0.0000, 0.0742, 0.1393, 0.2198, 0.3302, 0.5181, 1.0000},
        },
    };

    const Lserk *context = context_;
    number time_order = context->time_order;
    number num_stages = context->num_stages;
    const scalar *alpha = alpha_[time_order - 1][num_stages - 1];

    number num = eqns->mesh->cells.num_inner;
    number len = eqns->variables.len;
    number stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Update *primitive = eqns->variables.primitive;

    scalar(*variable)[stride] = eqns->variables.data;
    scalar(*derivative)[len] = arena_malloc(num, sizeof(*derivative));

    scalar full_timestep = courant * equations_timestep(eqns, variable, 0);
    scalar timestep = fmin(full_timestep, max_timestep);

    scalar(*variable0)[stride] = arena_memdup(variable, num, sizeof(*variable));
    for (number i = 0; i < num_stages; i++) {
        equations_derivative(eqns, variable, derivative, *time + (alpha[i] * timestep));
        for (number j = 0; j < num; j++) {
            for (number k = 0; k < len; k++) {
                variable[j][k] = variable0[j][k] + (derivative[j][k] * alpha[i + 1] * timestep);
            }
            primitive(variable[j], property);
        }
    }

    *time += timestep;
    equations_residual(eqns, derivative, residual_);

    arena_load(save);
    return full_timestep;
}
