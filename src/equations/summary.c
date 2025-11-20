#include <assert.h>
#include <string.h>

#include "equations.h"
#include "teal/utils.h"

void equations_summary(const Equations *eqns)
{
    assert(eqns);

    int len = 16;
    for (long i = 0; i < eqns->properties.num; i++) {
        len = lmax(len, strlen(eqns->properties.name[i]));
    }

    println("Equations summary");
    println("\t %-*s : %s", len, "name", eqns->name);

    for (long i = 0; i < eqns->properties.num; i++) {
        println("\t %-*s : %g", len, eqns->properties.name[i], eqns->properties.data[i]);
    }

    println("\t %-*s : %ld", len, "space order", eqns->space_order);
    if (*eqns->convective.name) {
        println("\t %-*s : %s", len, "convective flux", eqns->convective.name);
    }
    if (*eqns->viscous.name) {
        println("\t %-*s : %s", len, "viscous flux", eqns->viscous.name);
    }
    if (*eqns->limiter.name) {
        println("\t %-*s : %s", len, "limiter", eqns->limiter.name);
    }

    len = 0;
    for (long i = eqns->mesh->entities.num_inner; i < eqns->mesh->entities.num; i++) {
        len = lmax(len, strlen(eqns->mesh->entities.name[i]));
    }

    println("\t boundary conditions");
    for (long i = 0; i < eqns->boundary.num; i++) {
        println("\t\t %-*s : %s", len, eqns->boundary.entity[i], eqns->boundary.name[i]);
    }
    for (long i = eqns->mesh->entities.off_ghost; i < eqns->mesh->entities.num; i++) {
        println("\t\t %-*s : %s", len, eqns->mesh->entities.name[i], "periodic");
    }
}
