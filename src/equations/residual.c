#include <math.h>

#include "equations.h"
#include "teal/assert.h"
#include "teal/sync.h"
#include "teal/utils.h"

void equations_residual(const Equations *eqns, const void *derivative_, void *residual_)
{
    assert(eqns && derivative_);

    number num = eqns->mesh->cells.num_inner;
    scalar *volume = eqns->mesh->cells.volume;
    scalar sum_volume = eqns->mesh->cells.sum_volume;

    number len = eqns->variables.len;
    const scalar(*derivative)[len] = derivative_;
    scalar *residual = residual_;

    for (number i = 0; i < num; i++) {
        for (number j = 0; j < len; j++) {
            residual[j] += volume[i] * pow2(derivative[i][j]);
        }
    }
    for (number i = 0; i < len; i++) {
        residual[i] = sqrt(sync_fsum(residual[i])) / sum_volume;
    }
}
