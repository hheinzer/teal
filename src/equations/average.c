#include <assert.h>
#include <string.h>

#include "equations.h"
#include "sync.h"
#include "teal.h"

void equations_average(const Equations *eqns, const char *entity, const void *primitive_,
                       void *average_)
{
    assert(eqns && entity);

    double *volume = eqns->mesh->cells.volume;
    double sum_volume = eqns->mesh->cells.sum_volume;

    int num = eqns->mesh->entities.num_inner;
    String *name = eqns->mesh->entities.name;
    int *cell_off = eqns->mesh->entities.cell_off;

    int stride = eqns->primitive.stride;
    const double (*primitive)[stride] = primitive_;

    double *average = average_;
    memset(average, 0, stride * sizeof(*average));

    for (int i = 0; i < num; i++) {
        if (!strcmp(name[i], entity)) {
            for (int j = cell_off[i]; j < cell_off[i + 1]; j++) {
                for (int k = 0; k < stride; k++) {
                    average[k] += volume[j] * primitive[j][k];
                }
            }
            sync_sum(average, stride, MPI_DOUBLE);
            for (int k = 0; k < stride; k++) {
                average[k] /= sum_volume;
            }
            return;
        }
    }
    teal_error("invalid entity (%s)", entity);
}
