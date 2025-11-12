#include <math.h>

#include "equations.h"
#include "teal/arena.h"
#include "teal/assert.h"
#include "teal/utils.h"
#include "teal/vector.h"

scalar minmod(vector gradient, scalar variable, scalar minimum, scalar maximum, scalar parameter,
              const vector *offset, number beg, number end)
{
    unused(parameter);
    scalar psi = 1;
    for (number i = beg; i < end; i++) {
        scalar delta2 = vector_dot(gradient, offset[i]);
        if (isclose(delta2, 0)) {
            continue;
        }
        scalar delta1 = (delta2 > 0) ? (maximum - variable) : (minimum - variable);
        psi = fmin(psi, delta1 / delta2);
    }
    return psi;
}

scalar venkatakrishnan(vector gradient, scalar variable, scalar minimum, scalar maximum,
                       scalar parameter, const vector *offset, number beg, number end)
{
    scalar psi = 1;
    for (number i = beg; i < end; i++) {
        scalar delta2 = vector_dot(gradient, offset[i]);
        if (isclose(delta2, 0)) {
            continue;
        }
        scalar delta1 = (delta2 > 0) ? (maximum - variable) : (minimum - variable);
        scalar numerator = ((pow2(delta1) + parameter) * delta2) + (2 * pow2(delta2) * delta1);
        scalar denominator = pow2(delta1) + (2 * pow2(delta2)) + (delta1 * delta2) + parameter;
        psi = fmin(psi, numerator / denominator / delta2);
    }
    return psi;
}

scalar *venkatakrishnan_parameter(const Equations *eqns, scalar parameter)
{
    assert(eqns);
    number num = eqns->mesh->cells.num_inner;
    scalar *volume = eqns->mesh->cells.volume;
    scalar *eps2 = arena_malloc(num, sizeof(*eps2));
    for (number i = 0; i < num; i++) {
        eps2[i] = pow3(parameter * cbrt(volume[i]));
    }
    return eps2;
}

void equations_limiter(const Equations *eqns, const void *variable_, void *gradient_)
{
    assert(eqns && variable_ && gradient_);

    number num = eqns->mesh->cells.num_inner;
    number *cell_off = eqns->mesh->cells.cell.off;
    number *cell_idx = eqns->mesh->cells.cell.idx;
    vector *offset = eqns->mesh->cells.offset;

    number stride = eqns->variables.stride;
    scalar *parameter = eqns->limiter.parameter;
    Limiter *compute = eqns->limiter.compute;

    const scalar(*variable)[stride] = variable_;
    vector(*gradient)[stride] = gradient_;

    for (number i = 0; i < num; i++) {
        number beg = cell_off[i];
        number end = cell_off[i + 1];
        for (number j = 0; j < stride; j++) {
            scalar minimum = variable[i][j];
            scalar maximum = variable[i][j];
            for (number k = beg; k < end; k++) {
                minimum = fmin(minimum, variable[cell_idx[k]][j]);
                maximum = fmax(maximum, variable[cell_idx[k]][j]);
            }
            scalar psi = compute(gradient[i][j], variable[i][j], minimum, maximum,
                                 parameter ? parameter[i] : 0, offset, beg, end);
            vector_scale(&gradient[i][j], psi);
        }
    }
}
