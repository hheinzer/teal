#include <stdio.h>
#include <string.h>

#include "equations.h"
#include "teal/memory.h"
#include "teal/sync.h"

static long *compute_output_dimensions(String *name, long n_names);

Equations *equations_create(const Mesh *mesh, const char *name)
{
    memory_sum_setzero();

    Equations *eqns = memory_calloc(1, sizeof(*eqns));
    eqns->mesh = mesh;
    strcpy(eqns->name, name);

    const long n_inner_cells = mesh->n_inner_cells;
    eqns->timestep.value = memory_calloc(n_inner_cells, sizeof(*eqns->timestep.value));

    eqns->limiter.eps2 = memory_calloc(n_inner_cells, sizeof(*eqns->limiter.eps2));

    const long n_entities = mesh->n_entities;
    eqns->bc.name = memory_calloc(n_entities, sizeof(*eqns->bc.name));
    eqns->bc.state = memory_calloc(n_entities, sizeof(*eqns->bc.state));
    eqns->bc.apply = memory_calloc(n_entities, sizeof(*eqns->bc.apply));
    eqns->bc.compute = memory_calloc(n_entities, sizeof(*eqns->bc.compute));

    return eqns;
}

void equations_create_vars(Equations *eqns, const String *name, long n_cons, long n_vars)
{
    const long n_cells = eqns->mesh->n_cells;
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    eqns->n_cons = n_cons;
    eqns->n_vars = n_vars;
    eqns->vars.name = memory_duplicate(name, n_vars, sizeof(*name));
    eqns->vars.u = memory_calloc(n_cells * n_vars, sizeof(*eqns->vars.u));
    eqns->vars.dudt = memory_calloc(n_inner_cells * n_vars, sizeof(*eqns->vars.dudt));
    eqns->vars.dudx = memory_calloc(n_cells * n_vars, sizeof(*eqns->vars.dudx));
    eqns->vars.dim = compute_output_dimensions(eqns->vars.name, n_vars);

    const long n_send = eqns->mesh->sync.i_send[sync.size];
    eqns->sync.buf = memory_calloc(n_send * n_vars * N_DIMS, sizeof(*eqns->sync.buf));
    eqns->sync.recv = memory_calloc(sync.size, sizeof(*eqns->sync.recv));  // NOLINT
    eqns->sync.send = memory_calloc(sync.size, sizeof(*eqns->sync.send));  // NOLINT
}

void equations_create_scalar(Equations *eqns, const String *name, long n_scalars)
{
    eqns->n_scalars = n_scalars;
    eqns->scalar.name = memory_duplicate(name, n_scalars, sizeof(*name));
    eqns->scalar.value = memory_calloc(n_scalars, sizeof(*eqns->scalar.value));
}

void equations_create_user(Equations *eqns, const String *name, long n_user, Compute *compute)
{
    eqns->n_user = n_user;
    eqns->user.name = memory_duplicate(name, n_user, sizeof(*name));
    eqns->user.compute = compute;
    eqns->user.dim = compute_output_dimensions(eqns->user.name, n_user);
}

void equations_create_exact(Equations *eqns, Compute *compute)
{
    const long n_user = eqns->n_vars;
    smart String *name = memory_calloc(n_user, sizeof(*name));
    for (long i = 0; i < n_user; ++i) {
        sprintf(name[i], "exact %s", eqns->vars.name[i]);
        const long len = strlen(eqns->vars.name[i]);
        if (eqns->vars.name[i][len + 1]) {
            strcat(name[i], "-");
            strcat(name[i], &eqns->vars.name[i][len + 1]);
        }
    }
    equations_create_user(eqns, name, n_user, compute);
}

static long *compute_output_dimensions(String *name, long n_names)
{
    long *dim = memory_calloc(n_names, sizeof(*dim));
    for (long n = 0, i = 0; i < n_names; ++i) {
        dim[n] += 1;
        char *dash = strrchr(name[i], '-');
        if (dash) *dash = 0;
        if (i == n_names - 1 || strncmp(name[i + 1], name[i], strlen(name[i]))) n += 1;
    }
    return dim;
}
