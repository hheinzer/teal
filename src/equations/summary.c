#include <assert.h>

#include "equations.h"
#include "teal.h"

void equations_summary(const Equations *eqns)
{
    assert(eqns);

    int num_inner = eqns->mesh->entities.num_inner;
    String *entity = &eqns->mesh->entities.name[num_inner];

    int num_boundaries = eqns->boundary.num;
    String *boundary = eqns->boundary.name;

    int num_properties = eqns->properties.num;
    String *name = eqns->properties.name;
    double *property = eqns->properties.data;

    teal_print("Equations summary");
    teal_print("\t name            : %s", eqns->name);
    teal_print("\t space order     : %d", eqns->space_order);
    teal_print("\t convective flux : %s", eqns->convective.name);
    teal_print("\t viscous flux    : %s", eqns->viscous.name);
    teal_print("\t limiter         : %s", eqns->limiter.name);

    teal_print("\t Boundary conditions");
    for (int i = 0; i < num_boundaries; i++) {
        teal_print("\t\t %-20s : %s", entity[i], boundary[i]);
    }

    teal_print("\t Properties");
    for (int i = 0; i < num_properties; i++) {
        teal_print("\t\t %-20s : %g", name[i], property[i]);
    }
}
