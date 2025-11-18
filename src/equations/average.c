#include <string.h>

#include "equations.h"
#include "teal/arena.h"
#include "teal/assert.h"
#include "teal/sync.h"
#include "teal/utils.h"

void *equations_average(const Equations *eqns, const char *entity, const void *variable_)
{
    assert(eqns && entity && variable_);

    scalar *volume = eqns->mesh->cells.volume;
    scalar sum_volume = eqns->mesh->cells.sum_volume;

    int num = eqns->mesh->entities.num;
    Name *name = eqns->mesh->entities.name;
    int *cell_off = eqns->mesh->entities.cell_off;

    int stride = eqns->variables.stride;
    const scalar(*variable)[stride] = variable_;

    for (int i = 0; i < num; i++) {
        if (!strcmp(name[i], entity)) {
            scalar *average = arena_calloc(stride, sizeof(*average));
            for (int j = cell_off[i]; j < cell_off[i + 1]; j++) {
                for (int k = 0; k < stride; k++) {
                    average[k] += volume[j] * variable[j][k];
                }
            }
            for (int j = 0; j < stride; j++) {
                average[j] = sync_fsum(average[j]) / sum_volume;
            }
            return average;
        }
    }
    error("invalid entity -- '%s'", entity);
}
