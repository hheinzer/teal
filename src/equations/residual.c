#include <assert.h>
#include <math.h>
#include <string.h>

#include "equations.h"
#include "teal/sync.h"
#include "teal/utils.h"

void equations_residual(const Equations *eqns, const void *derivative_, void *residual_)
{
    assert(eqns && derivative_);

    long num = eqns->mesh->cells.num_inner;
    scalar *volume = eqns->mesh->cells.volume;
    scalar sum_volume = eqns->mesh->cells.sum_volume;

    long len = eqns->variables.len;
    const scalar(*derivative)[len] = derivative_;
    scalar *residual = residual_;

    memset(residual, 0, len * sizeof(*residual));
    for (long i = 0; i < num; i++) {
        for (long j = 0; j < len; j++) {
            residual[j] += volume[i] * sq(derivative[i][j]);
        }
    }
    for (long i = 0; i < len; i++) {
        residual[i] = sqrt(sync_fsum(residual[i]) / sum_volume);
    }
}
