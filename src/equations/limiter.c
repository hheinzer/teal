#include <assert.h>
#include <math.h>

#include "equations.h"
#include "teal/arena.h"
#include "teal/utils.h"
#include "teal/vector.h"

void vanleer(vector *gradient, scalar variable, scalar minimum, scalar maximum, scalar parameter,
             const vector *offset, long num)
{
    (void)parameter;
    scalar psi = 2;
    for (long i = 0; i < num; i++) {
        scalar delta2 = vector_dot(*gradient, offset[i]);
        if (delta2 == 0) {
            continue;
        }
        scalar delta1 = (delta2 > 0) ? (maximum - variable) : (minimum - variable);
        scalar ratio = delta1 / delta2;
        if (ratio <= 0) {
            psi = 0;
            break;
        }
        scalar phi = (2 * ratio) / (1 + ratio);
        psi = fmin(psi, phi);
    }
    vector_scale(gradient, psi);
}

void vanalbada1(vector *gradient, scalar variable, scalar minimum, scalar maximum, scalar parameter,
                const vector *offset, long num)
{
    (void)parameter;
    scalar psi = 1;
    for (long i = 0; i < num; i++) {
        scalar delta2 = vector_dot(*gradient, offset[i]);
        if (delta2 == 0) {
            continue;
        }
        scalar delta1 = (delta2 > 0) ? (maximum - variable) : (minimum - variable);
        scalar ratio = delta1 / delta2;
        if (ratio <= 0) {
            psi = 0;
            break;
        }
        scalar phi = (sq(ratio) + ratio) / (sq(ratio) + 1);
        psi = fmin(psi, phi);
    }
    vector_scale(gradient, psi);
}

void vanalbada2(vector *gradient, scalar variable, scalar minimum, scalar maximum, scalar parameter,
                const vector *offset, long num)
{
    (void)parameter;
    scalar psi = 1;
    for (long i = 0; i < num; i++) {
        scalar delta2 = vector_dot(*gradient, offset[i]);
        if (delta2 == 0) {
            continue;
        }
        scalar delta1 = (delta2 > 0) ? (maximum - variable) : (minimum - variable);
        scalar ratio = delta1 / delta2;
        if (ratio <= 0) {
            psi = 0;
            break;
        }
        scalar phi = (2 * ratio) / (sq(ratio) + 1);
        psi = fmin(psi, phi);
    }
    vector_scale(gradient, psi);
}

void mc(vector *gradient, scalar variable, scalar minimum, scalar maximum, scalar parameter,
        const vector *offset, long num)
{
    (void)parameter;
    scalar psi = 2;
    for (long i = 0; i < num; i++) {
        scalar delta2 = vector_dot(*gradient, offset[i]);
        if (delta2 == 0) {
            continue;
        }
        scalar delta1 = (delta2 > 0) ? (maximum - variable) : (minimum - variable);
        scalar ratio = delta1 / delta2;
        if (ratio <= 0) {
            psi = 0;
            break;
        }
        scalar phi = fmin(fmin(2 * ratio, (1 + ratio) / 2), 2);
        psi = fmin(psi, phi);
    }
    vector_scale(gradient, psi);
}

void koren(vector *gradient, scalar variable, scalar minimum, scalar maximum, scalar parameter,
           const vector *offset, long num)
{
    (void)parameter;
    scalar psi = 2;
    for (long i = 0; i < num; i++) {
        scalar delta2 = vector_dot(*gradient, offset[i]);
        if (delta2 == 0) {
            continue;
        }
        scalar delta1 = (delta2 > 0) ? (maximum - variable) : (minimum - variable);
        scalar ratio = delta1 / delta2;
        if (ratio <= 0) {
            psi = 0;
            break;
        }
        scalar phi = fmin(fmin(2 * ratio, (1 + 2 * ratio) / 3), 2);
        psi = fmin(psi, phi);
    }
    vector_scale(gradient, psi);
}

void minmod(vector *gradient, scalar variable, scalar minimum, scalar maximum, scalar parameter,
            const vector *offset, long num)
{
    (void)parameter;
    scalar psi = 1;
    for (long i = 0; i < num; i++) {
        scalar delta2 = vector_dot(*gradient, offset[i]);
        if (delta2 == 0) {
            continue;
        }
        scalar delta1 = (delta2 > 0) ? (maximum - variable) : (minimum - variable);
        scalar ratio = delta1 / delta2;
        if (ratio <= 0) {
            psi = 0;
            break;
        }
        scalar phi = fmin(1, ratio);
        psi = fmin(psi, phi);
    }
    vector_scale(gradient, psi);
}

void superbee(vector *gradient, scalar variable, scalar minimum, scalar maximum, scalar parameter,
              const vector *offset, long num)
{
    (void)parameter;
    scalar psi = 1;
    for (long i = 0; i < num; i++) {
        scalar delta2 = vector_dot(*gradient, offset[i]);
        if (delta2 == 0) {
            continue;
        }
        scalar delta1 = (delta2 > 0) ? (maximum - variable) : (minimum - variable);
        scalar ratio = delta1 / delta2;
        if (ratio <= 0) {
            psi = 0;
            break;
        }
        scalar phi = fmax(fmin(2 * ratio, 1), fmin(ratio, 2));
        psi = fmin(psi, phi);
    }
    vector_scale(gradient, psi);
}

void venkatakrishnan(vector *gradient, scalar variable, scalar minimum, scalar maximum,
                     scalar parameter, const vector *offset, long num)
{
    scalar psi = 1;
    for (long i = 0; i < num; i++) {
        scalar delta2 = vector_dot(*gradient, offset[i]);
        if (delta2 == 0) {
            continue;
        }
        scalar delta1 = (delta2 > 0) ? (maximum - variable) : (minimum - variable);
        scalar ratio = delta1 / delta2;
        if (ratio <= 0) {
            psi = 0;
            break;
        }
        scalar numerator = ((sq(delta1) + parameter) * delta2) + (2 * sq(delta2) * delta1);
        scalar denominator = sq(delta1) + (2 * sq(delta2)) + (delta1 * delta2) + parameter;
        psi = fmin(psi, numerator / denominator / delta2);
    }
    vector_scale(gradient, psi);
}

scalar *venkatakrishnan_parameter(const Equations *eqns, scalar parameter)
{
    assert(eqns);
    long num_inner = eqns->mesh->cells.num_inner;
    scalar *volume = eqns->mesh->cells.volume;
    scalar *eps2 = arena_malloc(num_inner, sizeof(*eps2));
    for (long i = 0; i < num_inner; i++) {
        eps2[i] = cb(parameter * cbrt(volume[i]));
    }
    return eps2;
}

void equations_limiter(const Equations *eqns, const void *variable_, void *gradient_)
{
    assert(eqns && variable_ && gradient_);

    long num_inner = eqns->mesh->cells.num_inner;
    long *cell_off = eqns->mesh->cells.cell.off;
    long *cell_idx = eqns->mesh->cells.cell.idx;
    vector *offset = eqns->mesh->cells.offset;

    long stride = eqns->variables.stride;
    scalar *parameter = eqns->limiter.parameter;
    Limiter *compute = eqns->limiter.compute;

    const scalar(*variable)[stride] = variable_;
    vector(*gradient)[stride] = gradient_;

    for (long i = 0; i < num_inner; i++) {
        long beg = cell_off[i];
        long end = cell_off[i + 1];
        for (long j = 0; j < stride; j++) {
            scalar minimum = variable[i][j];
            scalar maximum = variable[i][j];
            for (long k = beg; k < end; k++) {
                minimum = fmin(minimum, variable[cell_idx[k]][j]);
                maximum = fmax(maximum, variable[cell_idx[k]][j]);
            }
            compute(&gradient[i][j], variable[i][j], minimum, maximum, parameter ? parameter[i] : 0,
                    &offset[beg], end - beg);
        }
    }
}
