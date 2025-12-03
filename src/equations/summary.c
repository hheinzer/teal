#include <assert.h>
#include <string.h>

#include "equations.h"
#include "teal/utils.h"

void equations_summary(const Equations *eqns)
{
    assert(eqns);

    long num_entities = eqns->mesh->entities.num;
    long num_inner = eqns->mesh->entities.num_inner;
    long off_ghost = eqns->mesh->entities.off_ghost;
    Name *entity = eqns->mesh->entities.name;

    long num = eqns->boundary.num;
    Name *name = eqns->boundary.name;

    long width = 0;
    for (long i = num_inner; i < num_entities; i++) {
        width = lmax(width, strlen(entity[i]));
    }

    println("Equations summary");
    println("\t name            : %s", eqns->name);
    println("\t space order     : %ld", eqns->space_order);
    println("\t convective flux : %s", eqns->convective.name);
    println("\t viscous flux    : %s", eqns->viscous.name);
    println("\t limiter         : %s", eqns->limiter.name);

    println("\t Boundary conditions");
    for (long i = 0; i < num; i++) {
        println("\t\t %-*s : %s", width, entity[num_inner + i], name[i]);
    }
    for (long i = off_ghost; i < num_entities; i++) {
        println("\t\t %-*s : periodic", width, entity[i]);
    }
}
