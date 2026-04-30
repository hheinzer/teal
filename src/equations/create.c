#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "equations.h"
#include "teal.h"
#include "utils.h"

Equations *equations_create(const Mesh *mesh, const char *name, Timestep *timestep,
                            ConvectiveSelect *convective, ViscousSelect *viscous,
                            BoundarySelect *boundary)
{
    assert(mesh && name && timestep);

    Equations *eqns = teal_calloc(1, sizeof(*eqns));
    eqns->mesh = mesh;
    strcpy(eqns->name, name);
    eqns->timestep.compute = timestep;
    eqns->convective.select = convective;
    eqns->viscous.select = viscous;
    eqns->boundary.select = boundary;

    int num = mesh->entities.off_boundary - mesh->entities.num_inner;
    eqns->boundary.num = num;
    eqns->boundary.name = teal_calloc(num, sizeof(*eqns->boundary.name));
    eqns->boundary.context = teal_calloc(num, sizeof(*eqns->boundary.context));
    eqns->boundary.compute = teal_calloc(num, sizeof(*eqns->boundary.compute));
    eqns->boundary.condition = teal_calloc(num, sizeof(*eqns->boundary.condition));

    return eqns;
}

static void create_variables(EquationsVariables *variables, const Equations *eqns,
                             const int *dimension, int num)
{
    int stride = 0;
    for (int i = 0; i < num; i++) {
        assert(dimension[i] > 0);
        stride += dimension[i];
    }
    assert(stride > 0);

    variables->num = num;
    variables->stride = stride;

    variables->dimension = teal_calloc(num, sizeof(*dimension));
    copy(variables->dimension, dimension, num, sizeof(*dimension));

    variables->data = teal_calloc(eqns->mesh->cells.num, sizeof(double[stride]));
}

void equations_create_primitive(Equations *eqns, const char **name, const int *dimension,
                                Convert *convert, int num)
{
    assert(eqns && name && dimension && num > 0);
    eqns->primitive.name = teal_calloc(num, sizeof(*eqns->primitive.name));
    for (int i = 0; i < num; i++) {
        strcpy(eqns->primitive.name[i], name[i]);
    }
    create_variables(&eqns->primitive, eqns, dimension, num);
    eqns->primitive.func.convert = convert;
}

void equations_create_conserved(Equations *eqns, const char **name, const int *dimension,
                                Convert *convert, int num)
{
    assert(eqns && name && dimension && num > 0);
    eqns->conserved.name = teal_calloc(num, sizeof(*eqns->conserved.name));
    for (int i = 0; i < num; i++) {
        strcpy(eqns->conserved.name[i], name[i]);
    }
    create_variables(&eqns->conserved, eqns, dimension, num);
    eqns->conserved.func.convert = convert;
}

void equations_create_reference(Equations *eqns, Compute *compute)
{
    assert(eqns && compute);
    int num = eqns->primitive.num;
    eqns->reference.name = teal_calloc(num, sizeof(*eqns->reference.name));
    for (int i = 0; i < num; i++) {
        sprintf(eqns->reference.name[i], "%s-ref", eqns->primitive.name[i]);
    }
    create_variables(&eqns->reference, eqns, eqns->primitive.dimension, num);
    eqns->reference.func.compute = compute;
}

void equations_create_properties(Equations *eqns, const char **name, const double *property,
                                 int num)
{
    assert(eqns && name && property && num > 0);
    eqns->properties.num = num;
    eqns->properties.name = teal_calloc(num, sizeof(*eqns->properties.name));
    for (int i = 0; i < num; i++) {
        strcpy(eqns->properties.name[i], name[i]);
    }
    eqns->properties.data = teal_calloc(num, sizeof(*property));
    copy(eqns->properties.data, property, num, sizeof(*property));
}
