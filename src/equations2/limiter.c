#include "limiter.h"

#include <assert.h>
#include <math.h>

#include "teal2.h"

double *venkatakrishnan_parameter2(const Equations *eqns, double coefficient)
{
    assert(eqns);

    int num_inner = eqns->mesh->cells.num_inner;
    double *volume = eqns->mesh->cells.volume;

    double *eps2 = teal2_calloc(num_inner, sizeof(*eps2));
    for (int i = 0; i < num_inner; i++) {
        eps2[i] = cb(coefficient * cbrt(volume[i]));
    }

    return eps2;
}

static void barth_jespersen(const Equations *eqns, const void *primitive_, void *gradient_)
{
    int num_inner = eqns->mesh->cells.num_inner;
    int *cell_off = eqns->mesh->cells.cell.off;
    int *cell_idx = eqns->mesh->cells.cell.idx;
    Vector *offset = eqns->mesh->cells.offset;

    int stride = eqns->primitive.stride;
    const double (*primitive)[stride] = primitive_;
    Vector(*gradient)[stride] = gradient_;

    double min[stride];
    double max[stride];
    double psi[stride];

    for (int i = 0; i < num_inner; i++) {
        for (int k = 0; k < stride; k++) {
            min[k] = max[k] = primitive[i][k];
        }
        for (int j = cell_off[i]; j < cell_off[i + 1]; j++) {
            for (int k = 0; k < stride; k++) {
                min[k] = fmin(min[k], primitive[cell_idx[j]][k]);
                max[k] = fmax(max[k], primitive[cell_idx[j]][k]);
            }
        }
        for (int k = 0; k < stride; k++) {
            psi[k] = 1;
        }
        for (int j = cell_off[i]; j < cell_off[i + 1]; j++) {
            for (int k = 0; k < stride; k++) {
                double delta2 = vector2_dot(gradient[i][k], offset[j]);
                if (isclose(delta2, 0)) {
                    continue;
                }
                double delta1 = ((delta2 > 0) ? max[k] : min[k]) - primitive[i][k];
                psi[k] = fmin(psi[k], fmin(1, delta1 / delta2));
            }
        }
        for (int k = 0; k < stride; k++) {
            vector2_imul(&gradient[i][k], psi[k]);
        }
    }
}

static void venkatakrishnan(const Equations *eqns, const void *primitive_, void *gradient_)
{
    int num_inner = eqns->mesh->cells.num_inner;
    int *cell_off = eqns->mesh->cells.cell.off;
    int *cell_idx = eqns->mesh->cells.cell.idx;
    Vector *offset = eqns->mesh->cells.offset;

    double *eps2 = eqns->limiter.parameter;

    int stride = eqns->primitive.stride;
    const double (*primitive)[stride] = primitive_;
    Vector(*gradient)[stride] = gradient_;

    double min[stride];
    double max[stride];
    double psi[stride];

    for (int i = 0; i < num_inner; i++) {
        for (int k = 0; k < stride; k++) {
            min[k] = max[k] = primitive[i][k];
        }
        for (int j = cell_off[i]; j < cell_off[i + 1]; j++) {
            for (int k = 0; k < stride; k++) {
                min[k] = fmin(min[k], primitive[cell_idx[j]][k]);
                max[k] = fmax(max[k], primitive[cell_idx[j]][k]);
            }
        }
        for (int k = 0; k < stride; k++) {
            psi[k] = 1;
        }
        for (int j = cell_off[i]; j < cell_off[i + 1]; j++) {
            for (int k = 0; k < stride; k++) {
                double delta2 = vector2_dot(gradient[i][k], offset[j]);
                if (isclose(delta2, 0)) {
                    continue;
                }
                double delta1 = ((delta2 > 0) ? max[k] : min[k]) - primitive[i][k];
                double top = ((sq(delta1) + eps2[i]) * delta2) + (2 * sq(delta2) * delta1);
                double bot = sq(delta1) + (2 * sq(delta2)) + (delta1 * delta2) + eps2[i];
                psi[k] = fmin(psi[k], top / bot / delta2);
            }
        }
        for (int k = 0; k < stride; k++) {
            vector2_imul(&gradient[i][k], psi[k]);
        }
    }
}

void equations2_limiter(const Equations *eqns, const void *primitive, void *gradient)
{
    assert(eqns && primitive && gradient);
    switch (eqns->limiter.kind) {
        case BARTH_JESPERSEN: barth_jespersen(eqns, primitive, gradient); break;
        case VENKATAKRISHNAN: venkatakrishnan(eqns, primitive, gradient); break;
        default: teal2_error("invalid limiter kind (%d)", eqns->limiter.kind);
    }
}
