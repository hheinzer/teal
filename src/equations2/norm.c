#include <assert.h>
#include <math.h>
#include <string.h>

#include "equations2.h"
#include "sync2.h"

void equations2_norm(const Equations *eqns, void *norm_, double time)
{
    assert(eqns && norm_);

    int num_inner = eqns->mesh->cells.num_inner;
    double *volume = eqns->mesh->cells.volume;
    Vector *center = eqns->mesh->cells.center;
    double sum_volume = eqns->mesh->cells.sum_volume;

    int stride = eqns->primitive.stride;
    double (*primitive)[stride] = eqns->primitive.data;
    double (*reference)[stride] = eqns->reference.data;
    double *property = eqns->properties.data;
    Compute *compute = eqns->reference.fn.compute;

    double *norm = norm_;
    memset(norm, 0, stride * sizeof(*norm));

    for (int i = 0; i < num_inner; i++) {
        compute(reference[i], property, center[i], time, primitive[i]);
        for (int j = 0; j < stride; j++) {
            norm[j] += volume[i] * sq(reference[i][j] - primitive[i][j]);
        }
    }

    sync2_sum(norm, stride, MPI_DOUBLE);

    for (int i = 0; i < stride; i++) {
        norm[i] = sqrt(norm[i] / sum_volume);
    }
}
