#include <math.h>

#include "equations.h"
#include "teal/assert.h"

void equations_boundary(const Equations *eqns, void *variable_, scalar time)
{
    assert(eqns && variable_ && isfinite(time) && time >= 0);

    vector *center = eqns->mesh->cells.center;

    Adjacent *cell = eqns->mesh->faces.cell;
    matrix *basis = eqns->mesh->faces.basis;

    int stride = eqns->variables.stride;
    Update *conserved = eqns->variables.conserved;
    scalar *property = eqns->properties.data;

    int num = eqns->boundary.num;
    const int *face_off = eqns->boundary.face_off;
    const void **reference = eqns->boundary.reference;
    Compute **custom = eqns->boundary.custom;
    Boundary **condition = eqns->boundary.condition;

    scalar(*variable)[stride] = variable_;

    for (int i = 0; i < num; i++) {
        if (custom[i]) {
            for (int j = face_off[i]; j < face_off[i + 1]; j++) {
                int left = cell[j].left;
                int right = cell[j].right;
                custom[i](variable[right], property, center[right], time, variable[left]);
                conserved(variable[right], property);
            }
        }
        else if (condition[i]) {
            for (int j = face_off[i]; j < face_off[i + 1]; j++) {
                int left = cell[j].left;
                int right = cell[j].right;
                condition[i](variable[right], variable[left], reference[i], property, &basis[j]);
                conserved(variable[right], property);
            }
        }
    }
}
