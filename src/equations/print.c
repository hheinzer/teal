#include "teal/print.h"

#include <stdio.h>

#include "equations.h"
#include "teal/memory.h"
#include "teal/option.h"
#include "teal/sync.h"
#include "teal/utils.h"

void equations_print(const Equations *eqns)
{
    if (option.quiet) return;

    const long n_entities = eqns->mesh->n_entities;
    const alias(entity, eqns->mesh->entity.name);
    const alias(j_cell, eqns->mesh->entity.j_cell);
    long n_entity_cells[n_entities];
    for (long e = 0; e < n_entities; ++e) n_entity_cells[e] = sync_sum(j_cell[e + 1] - j_cell[e]);

    const double size = sync_sum(memory_sum_get());

    if (sync.rank == 0) {
        printf("%s equations summary:\n", eqns->name);

        for (long s = 0; s < eqns->n_scalars; ++s)
            print_key(eqns->scalar.name[s], "%g", eqns->scalar.value[s]);

        if (eqns->conv.flux) print_key("convective flux function", "%s", eqns->conv.name);
        if (eqns->visc.flux) print_key("viscous flux function", "%s", eqns->visc.name);

        print_key("space order", "%ld", eqns->space_order);
        print_key("limiter", "%s (k = %g)", eqns->limiter.name, eqns->limiter.k);

        for (long e = 0; e < n_entities; ++e)
            if (*eqns->bc.name[e])
                print_key("boundary condition", "%s(%ld) -> %s", entity[e], n_entity_cells[e],
                          eqns->bc.name[e]);

        print_size(size);
    }
}
