#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "equations.h"
#include "teal/arena.h"

Equations *equations_create(const Mesh *mesh, const char *name)
{
    assert(mesh && name);

    Equations *eqns = arena_calloc(1, sizeof(*eqns));

    eqns->mesh = mesh;
    strcpy(eqns->name, name);

    long num = mesh->entities.off_ghost - mesh->entities.num_inner;
    eqns->boundary.num = num;
    eqns->boundary.name = arena_calloc(num, sizeof(*eqns->boundary.name));
    eqns->boundary.reference = arena_calloc(num, sizeof(*eqns->boundary.reference));
    eqns->boundary.custom = arena_calloc(num, sizeof(*eqns->boundary.custom));
    eqns->boundary.condition = arena_calloc(num, sizeof(*eqns->boundary.condition));

    return eqns;
}

void equations_create_variables(Equations *eqns, const long *dim, const char **name,
                                Update *conserved, Update *primitive, long num_conserved, long num)
{
    assert(eqns && dim && name && conserved && primitive &&
           (0 < num_conserved && num_conserved <= num));

    long len = 0;
    long stride = 0;
    for (long i = 0; i < num; i++) {
        if (i < num_conserved) {
            len += dim[i];
        }
        stride += dim[i];
    }

    long num_cells = eqns->mesh->cells.num;
    scalar(*variable)[stride] = arena_calloc(num_cells, sizeof(*variable));

    eqns->variables.num = num;
    eqns->variables.len = len;
    eqns->variables.stride = stride;
    eqns->variables.dim = arena_memdup(dim, num, sizeof(*dim));
    eqns->variables.name = arena_calloc(num, sizeof(*eqns->variables.name));
    for (long i = 0; i < num; i++) {
        strcpy(eqns->variables.name[i], name[i]);
    }
    eqns->variables.data = variable;
    eqns->variables.conserved = conserved;
    eqns->variables.primitive = primitive;
}

void equations_create_properties(Equations *eqns, const char **name, const scalar *property,
                                 long num)
{
    assert(eqns && name && property && num > 0);
    eqns->properties.num = num;
    eqns->properties.name = arena_calloc(num, sizeof(*eqns->properties.name));
    for (long i = 0; i < num; i++) {
        strcpy(eqns->properties.name[i], name[i]);
    }
    eqns->properties.data = arena_memdup(property, num, sizeof(*property));
}

void equations_create_user_variables(Equations *eqns, const long *dim, const char **name,
                                     Compute *compute, long num)
{
    assert(eqns && dim && name && compute && num > 0);

    long stride = 0;
    for (long i = 0; i < num; i++) {
        stride += dim[i];
    }

    eqns->user.num = num;
    eqns->user.stride = stride;
    eqns->user.dim = arena_memdup(dim, num, sizeof(*dim));
    eqns->user.name = arena_calloc(num, sizeof(*eqns->user.name));
    for (long i = 0; i < num; i++) {
        strcpy(eqns->user.name[i], name[i]);
    }
    eqns->user.compute = compute;
}

void equations_create_exact_solution(Equations *eqns, Compute *compute)
{
    assert(eqns && compute);
    long num = eqns->variables.num;
    eqns->user.num = num;
    eqns->user.stride = eqns->variables.stride;
    eqns->user.dim = arena_memdup(eqns->variables.dim, num, sizeof(*eqns->variables.dim));
    eqns->user.name = arena_calloc(num, sizeof(*eqns->user.name));
    for (long i = 0; i < num; i++) {
        sprintf(eqns->user.name[i], "exact %s", eqns->variables.name[i]);
    }
    eqns->user.compute = compute;
    eqns->user.conserved = eqns->variables.conserved;
}
