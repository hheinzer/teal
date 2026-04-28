#include <assert.h>
#include <math.h>
#include <string.h>

#include "equations.h"
#include "sync.h"

void equations_residual(const Equations *eqns, const void *derivative_, void *residual_)
{
    assert(eqns && derivative_ && residual_);

    int num_inner = eqns->mesh->cells.num_inner;
    double *volume = eqns->mesh->cells.volume;
    double sum_volume = eqns->mesh->cells.sum_volume;

    int stride = eqns->conserved.stride;
    const double (*derivative)[stride] = derivative_;

    double *residual = residual_;
    memset(residual, 0, stride * sizeof(*residual));

    for (int i = 0; i < num_inner; i++) {
        for (int j = 0; j < stride; j++) {
            residual[j] += volume[i] * sq(derivative[i][j]);
        }
    }

    sync_sum(residual, stride, MPI_DOUBLE);

    for (int i = 0; i < stride; i++) {
        residual[i] = sqrt(residual[i] / sum_volume);
    }
}
