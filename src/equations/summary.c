#include "equations.h"
#include "teal/assert.h"
#include "teal/utils.h"

void equations_summary(const Equations *eqns)
{
    assert(eqns);

    println("Equations summary:");
    println("\t name:            %s", eqns->name);
    println("\t space order:     %td", eqns->space_order);
    println("\t convective flux: %s", optional(eqns->convective.name));
    println("\t viscous flux:    %s", optional(eqns->viscous.name));
    println("\t limiter:         %s", optional(eqns->limiter.name));

    println("\t boundary conditions:");
    for (number i = 0; i < eqns->boundary.num; i++) {
        println("\t\t %s: %s", eqns->boundary.entity[i], optional(eqns->boundary.name[i]));
    }
    for (number i = eqns->mesh->entities.off_ghost; i < eqns->mesh->entities.num; i++) {
        println("\t\t %s: %s", eqns->mesh->entities.name[i], "periodic");
    }
}
