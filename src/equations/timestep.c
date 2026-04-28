#include <assert.h>
#include <float.h>

#include "equations.h"
#include "sync.h"

double equations_timestep(const Equations *eqns, const void *primitive_)
{
    assert(eqns && primitive_);

    int num_inner = eqns->mesh->cells.num_inner;
    double *volume = eqns->mesh->cells.volume;
    Vector *projection = eqns->mesh->cells.projection;

    int stride = eqns->primitive.stride;
    const double (*primitive)[stride] = primitive_;
    double *property = eqns->properties.data;
    Timestep *compute = eqns->timestep.compute;

    double min_step = DBL_MAX;
    for (int i = 0; i < num_inner; i++) {
        double step = compute(primitive[i], property, volume[i], projection[i]);
        if (step < min_step) {
            min_step = step;
        }
    }

    sync_min(&min_step, 1, MPI_DOUBLE);
    return min_step;
}
