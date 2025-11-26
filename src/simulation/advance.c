#include "advance.h"

#include <math.h>

#include "gmres.h"
#include "teal/arena.h"
#include "teal/sync.h"

scalar euler(const Equations *eqns, scalar *time, void *residual_, scalar courant, scalar max_step,
             const void *ctx_)
{
    (void)ctx_;
    Arena save = arena_save();

    long num_inner = eqns->mesh->cells.num_inner;
    long len = eqns->variables.len;
    long stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Update *primitive = eqns->variables.primitive;

    scalar(*variable)[stride] = eqns->variables.data;
    scalar(*derivative)[len] = equations_derivative(eqns, variable, 0, *time);

    scalar step0 = courant * equations_time_step(eqns, variable, 0);
    scalar step = fmin(step0, max_step);

    for (long i = 0; i < num_inner; i++) {
        for (long j = 0; j < len; j++) {
            variable[i][j] += derivative[i][j] * step;
        }
        primitive(variable[i], property);
    }

    *time += step;
    equations_residual(eqns, derivative, residual_);

    arena_load(save);
    return step0;
}

scalar midpoint(const Equations *eqns, scalar *time, void *residual_, scalar courant,
                scalar max_step, const void *ctx_)
{
    (void)ctx_;
    Arena save = arena_save();

    long num_cells = eqns->mesh->cells.num;
    long num_inner = eqns->mesh->cells.num_inner;
    long len = eqns->variables.len;
    long stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Update *primitive = eqns->variables.primitive;

    scalar(*variable)[stride] = eqns->variables.data;
    scalar(*derivative)[len] = equations_derivative(eqns, variable, 0, *time);

    scalar step0 = courant * equations_time_step(eqns, variable, 0);
    scalar step = fmin(step0, max_step);

    scalar(*variable1)[stride] = arena_malloc(num_cells, sizeof(*variable1));
    for (long i = 0; i < num_inner; i++) {
        for (long j = 0; j < len; j++) {
            variable1[i][j] = variable[i][j] + (derivative[i][j] * step / 2);
        }
        primitive(variable1[i], property);
    }

    scalar(*derivative1)[len] = equations_derivative(eqns, variable1, 0, *time + (step / 2));

    for (long i = 0; i < num_inner; i++) {
        for (long j = 0; j < len; j++) {
            variable[i][j] += derivative1[i][j] * step;
        }
        primitive(variable[i], property);
    }

    *time += step;
    equations_residual(eqns, derivative, residual_);

    arena_load(save);
    return step0;
}

scalar heun(const Equations *eqns, scalar *time, void *residual_, scalar courant, scalar max_step,
            const void *ctx_)
{
    (void)ctx_;
    Arena save = arena_save();

    long num_cells = eqns->mesh->cells.num;
    long num_inner = eqns->mesh->cells.num_inner;
    long len = eqns->variables.len;
    long stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Update *primitive = eqns->variables.primitive;

    scalar(*variable)[stride] = eqns->variables.data;
    scalar(*derivative)[len] = equations_derivative(eqns, variable, 0, *time);

    scalar step0 = courant * equations_time_step(eqns, variable, 0);
    scalar step = fmin(step0, max_step);

    scalar(*variable1)[stride] = arena_malloc(num_cells, sizeof(*variable1));
    for (long i = 0; i < num_inner; i++) {
        for (long j = 0; j < len; j++) {
            variable1[i][j] = variable[i][j] + (derivative[i][j] * step);
        }
        primitive(variable1[i], property);
    }

    scalar(*derivative1)[len] = equations_derivative(eqns, variable1, 0, *time + step);

    for (long i = 0; i < num_inner; i++) {
        for (long j = 0; j < len; j++) {
            variable[i][j] += (derivative[i][j] + derivative1[i][j]) * step / 2;
        }
        primitive(variable[i], property);
    }

    *time += step;
    equations_residual(eqns, derivative, residual_);

    arena_load(save);
    return step0;
}

scalar ralston(const Equations *eqns, scalar *time, void *residual_, scalar courant,
               scalar max_step, const void *ctx_)
{
    (void)ctx_;
    Arena save = arena_save();

    long num_cells = eqns->mesh->cells.num;
    long num_inner = eqns->mesh->cells.num_inner;
    long len = eqns->variables.len;
    long stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Update *primitive = eqns->variables.primitive;

    scalar(*variable)[stride] = eqns->variables.data;
    scalar(*derivative)[len] = equations_derivative(eqns, variable, 0, *time);

    scalar step0 = courant * equations_time_step(eqns, variable, 0);
    scalar step = fmin(step0, max_step);

    scalar(*variable1)[stride] = arena_malloc(num_cells, sizeof(*variable1));
    for (long i = 0; i < num_inner; i++) {
        for (long j = 0; j < len; j++) {
            variable1[i][j] = variable[i][j] + (derivative[i][j] * 2 * step / 3);
        }
        primitive(variable1[i], property);
    }

    scalar(*derivative1)[len] = equations_derivative(eqns, variable1, 0, *time + (2 * step / 3));

    for (long i = 0; i < num_inner; i++) {
        for (long j = 0; j < len; j++) {
            variable[i][j] += step * (derivative[i][j] + 3 * derivative1[i][j]) / 4;
        }
        primitive(variable[i], property);
    }

    *time += step;
    equations_residual(eqns, derivative, residual_);

    arena_load(save);
    return step0;
}

scalar ssprk2(const Equations *eqns, scalar *time, void *residual_, scalar courant, scalar max_step,
              const void *ctx_)
{
    (void)ctx_;
    Arena save = arena_save();

    long num_cells = eqns->mesh->cells.num;
    long num_inner = eqns->mesh->cells.num_inner;
    long len = eqns->variables.len;
    long stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Update *primitive = eqns->variables.primitive;

    scalar(*variable)[stride] = eqns->variables.data;
    scalar(*derivative)[len] = equations_derivative(eqns, variable, 0, *time);

    scalar step0 = courant * equations_time_step(eqns, variable, 0);
    scalar step = fmin(step0, max_step);

    scalar(*variable1)[stride] = arena_malloc(num_cells, sizeof(*variable1));
    for (long i = 0; i < num_inner; i++) {
        for (long j = 0; j < len; j++) {
            variable1[i][j] = variable[i][j] + (derivative[i][j] * step);
        }
        primitive(variable1[i], property);
    }

    scalar(*derivative1)[len] = equations_derivative(eqns, variable1, 0, *time + step);

    for (long i = 0; i < num_inner; i++) {
        for (long j = 0; j < len; j++) {
            variable[i][j] = (variable[i][j] + variable1[i][j] + step * derivative1[i][j]) / 2;
        }
        primitive(variable[i], property);
    }

    *time += step;
    equations_residual(eqns, derivative, residual_);

    arena_load(save);
    return step0;
}

scalar ssprk3(const Equations *eqns, scalar *time, void *residual_, scalar courant, scalar max_step,
              const void *ctx_)
{
    (void)ctx_;
    Arena save = arena_save();

    long num_cells = eqns->mesh->cells.num;
    long num_inner = eqns->mesh->cells.num_inner;
    long len = eqns->variables.len;
    long stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Update *primitive = eqns->variables.primitive;

    scalar(*variable)[stride] = eqns->variables.data;
    scalar(*derivative)[len] = equations_derivative(eqns, variable, 0, *time);

    scalar step0 = courant * equations_time_step(eqns, variable, 0);
    scalar step = fmin(step0, max_step);

    scalar(*variable1)[stride] = arena_malloc(num_cells, sizeof(*variable1));
    for (long i = 0; i < num_inner; i++) {
        for (long j = 0; j < len; j++) {
            variable1[i][j] = variable[i][j] + (derivative[i][j] * step);
        }
        primitive(variable1[i], property);
    }

    scalar(*derivative1)[len] = equations_derivative(eqns, variable1, 0, *time + step);

    scalar(*variable2)[stride] = arena_malloc(num_cells, sizeof(*variable2));
    for (long i = 0; i < num_inner; i++) {
        for (long j = 0; j < len; j++) {
            variable2[i][j] = (3 * variable[i][j] + variable1[i][j] + step * derivative1[i][j]) / 4;
        }
        primitive(variable2[i], property);
    }

    scalar(*derivative2)[len] = equations_derivative(eqns, variable2, 0, *time + (step / 2));

    for (long i = 0; i < num_inner; i++) {
        for (long j = 0; j < len; j++) {
            variable[i][j] =
                (variable[i][j] + 2 * (variable2[i][j] + step * derivative2[i][j])) / 3;
        }
        primitive(variable[i], property);
    }

    *time += step;
    equations_residual(eqns, derivative, residual_);

    arena_load(save);
    return step0;
}

scalar rk3(const Equations *eqns, scalar *time, void *residual_, scalar courant, scalar max_step,
           const void *ctx_)
{
    (void)ctx_;
    Arena save = arena_save();

    long num_cells = eqns->mesh->cells.num;
    long num_inner = eqns->mesh->cells.num_inner;
    long len = eqns->variables.len;
    long stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Update *primitive = eqns->variables.primitive;

    scalar(*variable)[stride] = eqns->variables.data;
    scalar(*derivative)[len] = equations_derivative(eqns, variable, 0, *time);

    scalar step0 = courant * equations_time_step(eqns, variable, 0);
    scalar step = fmin(step0, max_step);

    scalar(*variable1)[stride] = arena_malloc(num_cells, sizeof(*variable1));
    for (long i = 0; i < num_inner; i++) {
        for (long j = 0; j < len; j++) {
            variable1[i][j] = variable[i][j] + (derivative[i][j] * step / 2);
        }
        primitive(variable1[i], property);
    }

    scalar(*derivative1)[len] = equations_derivative(eqns, variable1, 0, *time + (step / 2));

    scalar(*variable2)[stride] = arena_malloc(num_cells, sizeof(*variable2));
    for (long i = 0; i < num_inner; i++) {
        for (long j = 0; j < len; j++) {
            variable2[i][j] =
                variable[i][j] - (derivative[i][j] * step) + (derivative1[i][j] * 2 * step);
        }
        primitive(variable2[i], property);
    }

    scalar(*derivative2)[len] = equations_derivative(eqns, variable2, 0, *time + step);

    for (long i = 0; i < num_inner; i++) {
        for (long j = 0; j < len; j++) {
            variable[i][j] +=
                step * (derivative[i][j] + 4 * derivative1[i][j] + derivative2[i][j]) / 6;
        }
        primitive(variable[i], property);
    }

    *time += step;
    equations_residual(eqns, derivative, residual_);

    arena_load(save);
    return step0;
}

scalar rk4(const Equations *eqns, scalar *time, void *residual_, scalar courant, scalar max_step,
           const void *ctx_)
{
    (void)ctx_;
    Arena save = arena_save();

    long num_cells = eqns->mesh->cells.num;
    long num_inner = eqns->mesh->cells.num_inner;
    long len = eqns->variables.len;
    long stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Update *primitive = eqns->variables.primitive;

    scalar(*variable)[stride] = eqns->variables.data;
    scalar(*derivative)[len] = equations_derivative(eqns, variable, 0, *time);

    scalar step0 = courant * equations_time_step(eqns, variable, 0);
    scalar step = fmin(step0, max_step);

    scalar(*variable1)[stride] = arena_malloc(num_cells, sizeof(*variable1));
    for (long i = 0; i < num_inner; i++) {
        for (long j = 0; j < len; j++) {
            variable1[i][j] = variable[i][j] + (derivative[i][j] * step / 2);
        }
        primitive(variable1[i], property);
    }

    scalar(*derivative1)[len] = equations_derivative(eqns, variable1, 0, *time + (step / 2));

    scalar(*variable2)[stride] = arena_malloc(num_cells, sizeof(*variable2));
    for (long i = 0; i < num_inner; i++) {
        for (long j = 0; j < len; j++) {
            variable2[i][j] = variable[i][j] + (derivative1[i][j] * step / 2);
        }
        primitive(variable2[i], property);
    }

    scalar(*derivative2)[len] = equations_derivative(eqns, variable2, 0, *time + (step / 2));

    scalar(*variable3)[stride] = arena_malloc(num_cells, sizeof(*variable3));
    for (long i = 0; i < num_inner; i++) {
        for (long j = 0; j < len; j++) {
            variable3[i][j] = variable[i][j] + (derivative2[i][j] * step);
        }
        primitive(variable3[i], property);
    }

    scalar(*derivative3)[len] = equations_derivative(eqns, variable3, 0, *time + step);

    for (long i = 0; i < num_inner; i++) {
        for (long j = 0; j < len; j++) {
            variable[i][j] += step *
                              (derivative[i][j] + 2 * (derivative1[i][j] + derivative2[i][j]) +
                               derivative3[i][j]) /
                              6;
        }
        primitive(variable[i], property);
    }

    *time += step;
    equations_residual(eqns, derivative, residual_);

    arena_load(save);
    return step0;
}

scalar lserk(const Equations *eqns, scalar *time, void *residual_, scalar courant, scalar max_step,
             const void *ctx_)
{
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

    const RungeKutta *ctx = ctx_;
    long time_order = ctx->time_order;
    long num_stages = ctx->num_stages;
    const scalar *alpha = alpha_[time_order - 1][num_stages - 1];

    long num_inner = eqns->mesh->cells.num_inner;
    long len = eqns->variables.len;
    long stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Update *primitive = eqns->variables.primitive;

    scalar(*variable)[stride] = eqns->variables.data;
    scalar(*variable0)[stride] = arena_memdup(variable, num_inner, sizeof(*variable));

    scalar step0 = courant * equations_time_step(eqns, variable, 0);
    scalar step = fmin(step0, max_step);

    scalar(*derivative)[len] = arena_malloc(num_inner, sizeof(*derivative));
    for (long i = 0; i < num_stages; i++) {
        equations_derivative(eqns, variable, derivative, *time + (alpha[i] * step));
        for (long j = 0; j < num_inner; j++) {
            for (long k = 0; k < len; k++) {
                variable[j][k] = variable0[j][k] + (derivative[j][k] * alpha[i + 1] * step);
            }
            primitive(variable[j], property);
        }
    }

    *time += step;
    equations_residual(eqns, derivative, residual_);

    arena_load(save);
    return step0;
}

scalar implicit_euler(const Equations *eqns, scalar *time, void *residual_, scalar courant,
                      scalar max_step, const void *ctx_)
{
    Arena save = arena_save();

    const NewtonKrylov *ctx = ctx_;
    scalar tol = ctx->newton_tolerance;

    long num_inner = eqns->mesh->cells.num_inner;
    long len = eqns->variables.len;
    long stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Update *primitive = eqns->variables.primitive;

    scalar(*variable)[stride] = eqns->variables.data;
    scalar(*variable0)[stride] = arena_memdup(variable, num_inner, sizeof(*variable));

    scalar step0 = courant * equations_time_step(eqns, variable, 0);
    scalar step = fmin(step0, max_step);

    *time += step;
    scalar(*derivative)[len] = equations_derivative(eqns, variable, 0, *time);

    scalar(*residual)[len] = arena_malloc(num_inner, sizeof(*residual));
    for (long i = 0; i < num_inner; i++) {
        for (long j = 0; j < len; j++) {
            residual[i][j] = -derivative[i][j] * step;
        }
    }

    scalar norm = sync_norm(*residual, num_inner * len);

    static const scalar fd_relative = 1e-6;
    scalar fd_scale = fd_relative * fmax(1, sync_norm(*variable, num_inner * len)) / fmax(1, norm);

    scalar tol_norm = tol * norm;
    long max_iter = 128 * (num_inner + 1);

    scalar(*increment)[len] = arena_malloc(num_inner, sizeof(*increment));
    for (long iter = 0; iter < max_iter && norm > tol_norm; iter++) {
        gmres(eqns, variable, derivative, residual, increment, *time, step, norm, fd_scale, ctx);
        for (long i = 0; i < num_inner; i++) {
            for (long j = 0; j < len; j++) {
                variable[i][j] += increment[i][j];
            }
            primitive(variable[i], property);
        }

        equations_derivative(eqns, variable, derivative, *time);
        for (long i = 0; i < num_inner; i++) {
            for (long j = 0; j < len; j++) {
                residual[i][j] = variable[i][j] - variable0[i][j] - (derivative[i][j] * step);
            }
        }

        norm = sync_norm(*residual, num_inner * len);

        fd_scale = fd_relative * fmax(1, sync_norm(*variable, num_inner * len)) / fmax(1, norm);
    }

    equations_residual(eqns, derivative, residual_);

    arena_load(save);
    return step0;
}
