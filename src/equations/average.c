#include <assert.h>
#include <string.h>

#include "equations.h"
#include "teal/sync.h"
#include "teal/utils.h"

void equations_average(const Equations *eqns, const char *entity, void *average_)
{
    assert(eqns && entity);

    scalar *volume = eqns->mesh->cells.volume;
    scalar sum_volume = eqns->mesh->cells.sum_volume;

    long num = eqns->mesh->entities.num;
    Name *name = eqns->mesh->entities.name;
    long *cell_off = eqns->mesh->entities.cell_off;

    long stride = eqns->variables.stride;
    const scalar(*variable)[stride] = eqns->variables.data;

    scalar *average = average_;
    memset(average, 0, stride * sizeof(*average));

    for (long i = 0; i < num; i++) {
        if (!strcmp(name[i], entity)) {
            for (long j = cell_off[i]; j < cell_off[i + 1]; j++) {
                for (long k = 0; k < stride; k++) {
                    average[k] += volume[j] * variable[j][k];
                }
            }
            for (long j = 0; j < stride; j++) {
                average[j] = sync_fsum(average[j]) / sum_volume;
            }
            return;
        }
    }
    error("invalid entity (%s)", entity);
}
