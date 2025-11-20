#include <assert.h>
#include <math.h>

#include "equations.h"
#include "teal/arena.h"
#include "teal/utils.h"

scalar minmod(vector gradient, scalar variable, scalar minimum, scalar maximum, scalar parameter,
              const vector *offset, long beg, long end)
{
    (void)parameter;
    scalar psi = 1;
    for (long i = beg; i < end; i++) {
        scalar delta2 =
            (gradient.x * offset[i].x) + (gradient.y * offset[i].y) + (gradient.z * offset[i].z);
        if (delta2 == 0) {
            continue;
        }
        scalar delta1 = (delta2 > 0) ? (maximum - variable) : (minimum - variable);
        psi = fmin(psi, delta1 / delta2);
    }
    return psi;
}

scalar venkatakrishnan(vector gradient, scalar variable, scalar minimum, scalar maximum,
                       scalar parameter, const vector *offset, long beg, long end)
{
    scalar psi = 1;
    for (long i = beg; i < end; i++) {
        scalar delta2 =
            (gradient.x * offset[i].x) + (gradient.y * offset[i].y) + (gradient.z * offset[i].z);
        if (delta2 == 0) {
            continue;
        }
        scalar delta1 = (delta2 > 0) ? (maximum - variable) : (minimum - variable);
        scalar delta12 = sq(delta1);
        scalar delta22 = sq(delta2);
        scalar numerator = ((delta12 + parameter) * delta2) + (2 * delta22 * delta1);
        scalar denominator = delta12 + (2 * delta22) + (delta1 * delta2) + parameter;
        psi = fmin(psi, numerator / denominator / delta2);
    }
    return psi;
}

scalar *venkatakrishnan_parameter(const Equations *eqns, scalar parameter)
{
    assert(eqns);
    long num = eqns->mesh->cells.num_inner;
    scalar *volume = eqns->mesh->cells.volume;
    scalar *eps2 = arena_malloc(num, sizeof(*eps2));
    for (long i = 0; i < num; i++) {
        eps2[i] = cb(parameter * cbrt(volume[i]));
    }
    return eps2;
}

void equations_limiter(const Equations *eqns, const void *variable_, void *gradient_)
{
    assert(eqns && variable_ && gradient_);

    long num = eqns->mesh->cells.num_inner;
    long *cell_off = eqns->mesh->cells.cell.off;
    long *cell_idx = eqns->mesh->cells.cell.idx;
    vector *offset = eqns->mesh->cells.offset;

    long stride = eqns->variables.stride;
    scalar *parameter = eqns->limiter.parameter;
    Limiter *compute = eqns->limiter.compute;

    const scalar(*variable)[stride] = variable_;
    vector(*gradient)[stride] = gradient_;

    for (long i = 0; i < num; i++) {
        long beg = cell_off[i];
        long end = cell_off[i + 1];
        for (long j = 0; j < stride; j++) {
            scalar minimum = variable[i][j];
            scalar maximum = variable[i][j];
            for (long k = beg; k < end; k++) {
                minimum = fmin(minimum, variable[cell_idx[k]][j]);
                maximum = fmax(maximum, variable[cell_idx[k]][j]);
            }
            scalar psi = compute(gradient[i][j], variable[i][j], minimum, maximum,
                                 parameter ? parameter[i] : 0, offset, beg, end);
            gradient[i][j].x *= psi;
            gradient[i][j].y *= psi;
            gradient[i][j].z *= psi;
        }
    }
}
