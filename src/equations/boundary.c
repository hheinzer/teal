#include <assert.h>
#include <math.h>

#include "equations.h"

void equations_boundary(const Equations *eqns, void *variable_, scalar time)
{
    assert(eqns && variable_ && isfinite(time) && time >= 0);

    vector *center = eqns->mesh->cells.center;

    Adjacent *cell = eqns->mesh->faces.cell;
    Basis *basis = eqns->mesh->faces.basis;

    long num_inner = eqns->mesh->entities.num_inner;
    long *face_off = &eqns->mesh->entities.face_off[num_inner];

    long stride = eqns->variables.stride;
    Update *conserved = eqns->variables.conserved;
    scalar *property = eqns->properties.data;

    long num = eqns->boundary.num;
    const void **reference = eqns->boundary.reference;
    Compute **custom = eqns->boundary.custom;
    Boundary **condition = eqns->boundary.condition;

    scalar(*variable)[stride] = variable_;

    for (long i = 0; i < num; i++) {
        if (custom[i]) {
            for (long j = face_off[i]; j < face_off[i + 1]; j++) {
                long left = cell[j].left;
                long right = cell[j].right;
                custom[i](variable[right], property, center[right], time, variable[left]);
                conserved(variable[right], property);
            }
        }
        else if (condition[i]) {
            for (long j = face_off[i]; j < face_off[i + 1]; j++) {
                long left = cell[j].left;
                long right = cell[j].right;
                condition[i](variable[right], variable[left], reference[i], property, &basis[j]);
                conserved(variable[right], property);
            }
        }
    }
}
