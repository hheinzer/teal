#include <assert.h>
#include <string.h>

#include "equations.h"
#include "teal/utils.h"

static long max_length(const Name *name, long beg, long end)
{
    long max = 0;
    for (long i = beg; i < end; i++) {
        max = lmax(max, strlen(name[i]));
    }
    return max;
}

void equations_summary(const Equations *eqns)
{
    assert(eqns);

    long num_entities = eqns->mesh->entities.num;
    long num_inner = eqns->mesh->entities.num_inner;
    long off_ghost = eqns->mesh->entities.off_ghost;
    Name *entity = eqns->mesh->entities.name;

    long num_properties = eqns->properties.num;
    long num_boundaries = eqns->boundary.num;
    Name *property = eqns->properties.name;
    scalar *value = eqns->properties.data;
    Name *name = eqns->boundary.name;

    println("Equations summary");
    println("\t name            : %s", eqns->name);
    println("\t space order     : %ld", eqns->space_order);
    println("\t convective flux : %s", eqns->convective.name);
    println("\t viscous flux    : %s", eqns->viscous.name);
    println("\t limiter         : %s", eqns->limiter.name);

    println("\t Material properties");
    int width = max_length(property, 0, num_properties);
    for (long i = 0; i < num_properties; i++) {
        println("\t\t %-*s : %g", width, property[i], value[i]);
    }

    println("\t Boundary conditions");
    width = max_length(entity, num_inner, num_entities);
    for (long i = 0; i < num_boundaries; i++) {
        println("\t\t %-*s : %s", width, entity[num_inner + i], name[i]);
    }
    for (long i = off_ghost; i < num_entities; i++) {
        println("\t\t %-*s : periodic", width, entity[i]);
    }
}
