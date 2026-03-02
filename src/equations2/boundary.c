#include <assert.h>

#include "equations2.h"

void equations2_boundary(const Equations *eqns, void *primitive_, double time)
{
    assert(eqns && primitive_);

    Vector *center = eqns->mesh->cells.center;

    Pair *cell_idx = eqns->mesh->faces.cell_idx;
    Matrix *basis = eqns->mesh->faces.basis;

    int num_inner = eqns->mesh->entities.num_inner;
    int *face_off = &eqns->mesh->entities.face_off[num_inner];

    int stride = eqns->primitive.stride;
    double (*primitive)[stride] = primitive_;
    double *property = eqns->properties.data;

    int num = eqns->boundary.num;
    Compute **compute = eqns->boundary.compute;
    Boundary **condition = eqns->boundary.condition;
    const void **context = eqns->boundary.context;

    for (int i = 0; i < num; i++) {
        if (compute[i]) {
            for (int j = face_off[i]; j < face_off[i + 1]; j++) {
                int left = cell_idx[j].left;
                int right = cell_idx[j].right;
                compute[i](primitive[right], property, center[right], time, primitive[left]);
            }
        }
        else {
            assert(condition[i]);
            for (int j = face_off[i]; j < face_off[i + 1]; j++) {
                int left = cell_idx[j].left;
                int right = cell_idx[j].right;
                condition[i](primitive[right], primitive[left], property, basis[j], context[i]);
            }
        }
    }
}
