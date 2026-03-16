#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "equations2.h"
#include "teal2.h"
#include "utils2.h"

static void wrap_destroy(void *eqns)
{
    equations2_destroy(eqns);
}

Equations *equations2_create(const Mesh *mesh, const char *name, Timestep *timestep,
                             ConvectiveSelect *convective, ViscousSelect *viscous,
                             BoundarySelect *boundary, int space_order)
{
    assert(mesh && name && 1 <= space_order && space_order <= 2);

    Equations *eqns = teal2_calloc(1, sizeof(*eqns));

    eqns->mesh = mesh;
    strcpy(eqns->name, name);
    eqns->space_order = space_order;

    eqns->timestep.compute = timestep;
    eqns->boundary.select = boundary;
    eqns->convective.select = convective;
    eqns->viscous.select = viscous;

    int num = mesh->entities.off_boundary - mesh->entities.num_inner;
    eqns->boundary.num = num;
    eqns->boundary.name = teal2_calloc(num, sizeof(*eqns->boundary.name));
    eqns->boundary.context = teal2_calloc(num, sizeof(*eqns->boundary.context));
    eqns->boundary.compute = teal2_calloc(num, sizeof(*eqns->boundary.compute));
    eqns->boundary.condition = teal2_calloc(num, sizeof(*eqns->boundary.condition));

    return eqns;
}

static void create_variables(EquationsVariables *variables, const Equations *eqns,
                             const int *dimension, Compute *compute, int num)
{
    int stride = 0;
    for (int i = 0; i < num; i++) {
        assert(dimension[i] > 0);
        stride += dimension[i];
    }
    assert(stride > 0);

    variables->num = num;
    variables->stride = stride;

    variables->dimension = teal2_calloc(num, sizeof(*dimension));
    copy(variables->dimension, dimension, num, sizeof(*dimension));

    variables->data = teal2_calloc(eqns->mesh->cells.num, (int)sizeof(double[stride]));
    variables->compute = compute;
}

void equations2_create_primitive(Equations *eqns, const char **name, const int *dimension,
                                 Compute *compute, int num)
{
    assert(eqns && name && dimension && compute && num > 0);

    eqns->primitive.name = teal2_calloc(num, sizeof(*eqns->primitive.name));
    for (int i = 0; i < num; i++) {
        strcpy(eqns->primitive.name[i], name[i]);
    }

    create_variables(&eqns->primitive, eqns, dimension, compute, num);
}

void equations2_create_conserved(Equations *eqns, const char **name, const int *dimension,
                                 Compute *compute, int num)
{
    assert(eqns && name && dimension && compute && num > 0);

    eqns->conserved.name = teal2_calloc(num, sizeof(*eqns->conserved.name));
    for (int i = 0; i < num; i++) {
        strcpy(eqns->conserved.name[i], name[i]);
    }

    create_variables(&eqns->conserved, eqns, dimension, compute, num);
}

void equations2_create_reference(Equations *eqns, Compute *compute)
{
    assert(eqns && compute);

    int num = eqns->primitive.num;

    eqns->reference.name = teal2_calloc(num, sizeof(*eqns->reference.name));
    for (int i = 0; i < eqns->primitive.num; i++) {
        sprintf(eqns->reference.name[i], "%s-ref", eqns->primitive.name[i]);
    }

    create_variables(&eqns->reference, eqns, eqns->primitive.dimension, compute, num);
}

void equations2_create_properties(Equations *eqns, const char **name, const double *property,
                                  int num)
{
    assert(eqns && name && property && num > 0);

    eqns->properties.num = num;

    eqns->properties.name = teal2_calloc(num, sizeof(*eqns->properties.name));
    for (int i = 0; i < num; i++) {
        strcpy(eqns->properties.name[i], name[i]);
    }

    eqns->properties.data = teal2_calloc(num, sizeof(*property));
    copy(eqns->properties.data, property, num, sizeof(*property));
}
