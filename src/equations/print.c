#include <stdio.h>

#include "core/sync.h"
#include "core/utils.h"
#include "equations.h"
#include "teal.h"

static double compute_size(const Equations *eqns);

void equations_print(const Equations *eqns)
{
    const long n_entities = eqns->mesh->n_entities;
    const ALIAS(entity, eqns->mesh->entity.name);
    const ALIAS(j_cell, eqns->mesh->entity.j_cell);
    long n_entity_cells[n_entities];
    for (long e = 0; e < n_entities; ++e) n_entity_cells[e] = sync_sum(j_cell[e + 1] - j_cell[e]);

    double size = sync_sum(compute_size(eqns));
    const char mod = sizefmt(&size);

    if (teal.rank == 0) {
        printf("%s equations summary:\n", eqns->name);
        printf(" | " KEYFMT ": %ld\n", "space order", eqns->space_order);
        for (long s = 0; s < eqns->n_scalars; ++s)
            printf(" | " KEYFMT ": %g\n", eqns->scalar.name[s], eqns->scalar.value[s]);
        printf(" | " KEYFMT ": %s\n", "convective flux function", eqns->flux.name_conv);
        if (*eqns->limiter.name)
            printf(" | " KEYFMT ": %s (k = %g)\n", "limiter", eqns->limiter.name, eqns->limiter.k);
        for (long e = 0; e < n_entities; ++e)
            if (*eqns->bc.name[e])
                printf(" | " KEYFMT ": %s(%ld) -> %s\n", "boundary condition", entity[e],
                       n_entity_cells[e], eqns->bc.name[e]);
        printf(" | " KEYFMT ": %g %cB\n", "memory size", size, mod);
    }
}

static double compute_size(const Equations *eqns)
{
    double size = sizeof(*eqns);

    size += eqns->n_scalars * sizeof(*eqns->scalar.name);
    size += eqns->n_scalars * sizeof(*eqns->scalar.value);

    const long n_cells = eqns->mesh->n_cells;
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    size += eqns->n_vars * sizeof(*eqns->vars.name);
    size += n_cells * eqns->n_vars * sizeof(*eqns->vars.u);
    if (eqns->space_order == 2) size += n_cells * eqns->n_vars * N_DIMS * sizeof(*eqns->vars.dudx);
    size += n_inner_cells * eqns->n_vars * sizeof(*eqns->vars.dudt);
    size += eqns->n_vars * sizeof(*eqns->vars.dim);

    if (eqns->limiter.func) size += n_inner_cells * sizeof(*eqns->limiter.eps2);

    size += eqns->n_user * sizeof(*eqns->user.name);
    size += eqns->n_user * sizeof(*eqns->user.dim);

    const long n_entities = eqns->mesh->n_entities;
    size += n_entities * sizeof(*eqns->bc.name);
    size += n_entities * sizeof(*eqns->bc.state);
    size += n_entities * sizeof(*eqns->bc.apply);
    size += n_entities * sizeof(*eqns->bc.custom);

    const long n_send = eqns->mesh->sync.i_send[teal.size];
    size += n_send * eqns->n_vars * N_DIMS * sizeof(*eqns->sync.buf);
    size += teal.size * sizeof(*eqns->sync.recv);  // NOLINT
    size += teal.size * sizeof(*eqns->sync.send);  // NOLINT

    return size;
}
