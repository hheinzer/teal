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

    int num = mesh->entities.off_ghost - mesh->entities.num_inner;
    eqns->boundary.num = num;
    eqns->boundary.entity = &mesh->entities.name[mesh->entities.num_inner];
    eqns->boundary.cell_off = &mesh->entities.cell_off[mesh->entities.num_inner];
    eqns->boundary.face_off = &mesh->entities.face_off[mesh->entities.num_inner];
    eqns->boundary.name = arena_calloc(num, sizeof(*eqns->boundary.name));
    eqns->boundary.reference = arena_calloc(num, sizeof(*eqns->boundary.reference));
    eqns->boundary.custom = arena_calloc(num, sizeof(*eqns->boundary.custom));
    eqns->boundary.condition = arena_calloc(num, sizeof(*eqns->boundary.condition));

    return eqns;
}

void equations_create_variables(Equations *eqns, const char **name, const Type *type,
                                Update *conserved, Update *primitive, int num_conserved, int num)
{
    assert(eqns && name && type && 0 < num_conserved && num_conserved <= num);
    int len = 0;
    int stride = 0;
    for (int i = 0; i < num; i++) {
        if (i < num_conserved) {
            len += type[i];
        }
        stride += type[i];
    }
    eqns->variables.num = num;
    eqns->variables.len = len;
    eqns->variables.stride = stride;
    eqns->variables.name = arena_calloc(num, sizeof(*eqns->variables.name));
    for (int i = 0; i < num; i++) {
        strcpy(eqns->variables.name[i], name[i]);
    }
    eqns->variables.type = arena_memdup(type, num, sizeof(*type));
    eqns->variables.data = arena_calloc(eqns->mesh->cells.num * stride, sizeof(scalar));
    eqns->variables.conserved = conserved;
    eqns->variables.primitive = primitive;
}

void equations_create_properties(Equations *eqns, const char **name, const scalar *property,
                                 int num)
{
    assert(eqns && name && property && num > 0);
    eqns->properties.num = num;
    eqns->properties.name = arena_calloc(num, sizeof(*eqns->properties.name));
    for (int i = 0; i < num; i++) {
        strcpy(eqns->properties.name[i], name[i]);
    }
    eqns->properties.data = arena_memdup(property, num, sizeof(*property));
}

void equations_create_user_variables(Equations *eqns, const char **name, const Type *type,
                                     Compute *compute, int num)
{
    assert(eqns && name && type && compute && num > 0);
    int stride = 0;
    for (int i = 0; i < num; i++) {
        stride += type[i];
    }
    eqns->user.num = num;
    eqns->user.stride = stride;
    eqns->user.name = arena_calloc(num, sizeof(*eqns->user.name));
    for (int i = 0; i < num; i++) {
        strcpy(eqns->user.name[i], name[i]);
    }
    eqns->user.type = arena_memdup(type, num, sizeof(*type));
    eqns->user.compute = compute;
}

void equations_create_exact_solution(Equations *eqns, Compute *compute)
{
    assert(eqns && compute);
    int num = eqns->variables.num;
    eqns->user.num = num;
    eqns->user.stride = eqns->variables.stride;
    eqns->user.name = arena_calloc(num, sizeof(*eqns->user.name));
    for (int i = 0; i < num; i++) {
        sprintf(eqns->user.name[i], "exact %s", eqns->variables.name[i]);
    }
    eqns->user.type = arena_memdup(eqns->variables.type, num, sizeof(*eqns->variables.type));
    eqns->user.compute = compute;
    eqns->user.conserved = eqns->variables.conserved;
}
