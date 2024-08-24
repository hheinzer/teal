#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core/memory.h"
#include "core/utils.h"
#include "equations.h"
#include "limiter.h"
#include "teal.h"

static void compute_output_dimensions(char (*name)[NAMELEN], long **dim, long n_names);

void equations_create(Equations *eqns, long space_order)
{
    if (space_order < 1 || 2 < space_order) error("unsupported space order '%ld'", space_order);
    eqns->space_order = space_order;

    const long n_cons = eqns->n_cons;
    const long n_vars = eqns->n_vars;
    const long n_cells = eqns->mesh->n_cells;
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    eqns->vars.u = memory_calloc(n_cells * n_vars, sizeof(*eqns->vars.u));
    eqns->vars.dudt = memory_calloc(n_inner_cells * n_cons, sizeof(*eqns->vars.dudt));
    compute_output_dimensions(eqns->vars.name, &eqns->vars.dim, eqns->n_vars);

    if (eqns->space_order == 2) {
        eqns->vars.dudx = memory_calloc(n_cells * n_vars * N_DIMS, sizeof(*eqns->vars.dudx));
        equations_set_limiter(eqns, "minmod", 0);
    }

    const long n_entities = eqns->mesh->n_entities;
    eqns->bc.name = memory_calloc(n_entities, sizeof(*eqns->bc.name));
    eqns->bc.state = memory_calloc(n_entities, sizeof(*eqns->bc.state));
    eqns->bc.apply = memory_calloc(n_entities, sizeof(*eqns->bc.apply));
    eqns->bc.func = memory_calloc(n_entities, sizeof(*eqns->bc.func));

    const long n_send = eqns->mesh->sync.i_send[teal.size];
    eqns->sync.buf = memory_calloc(n_send * n_vars * N_DIMS, sizeof(*eqns->sync.buf));
    eqns->sync.recv = memory_calloc(teal.size, sizeof(*eqns->sync.recv));  // NOLINT
    eqns->sync.send = memory_calloc(teal.size, sizeof(*eqns->sync.send));  // NOLINT
}

void equations_create_user(Equations *eqns, const char **name, Function *func, long n_user)
{
    eqns->n_user = n_user;
    eqns->user.name = memory_calloc(n_user, sizeof(*eqns->user.name));
    for (long v = 0; v < n_user; ++v) strcpy(eqns->user.name[v], name[v]);
    eqns->user.func = func;
    compute_output_dimensions(eqns->user.name, &eqns->user.dim, eqns->n_user);
}

void equations_create_exact(Equations *eqns, Function *func)
{
    const long n_user = eqns->n_vars;
    smart char **name = memory_calloc(n_user, sizeof(*name));
    for (long v = 0; v < n_user; ++v) {
        char buf[NAMELEN];
        sprintf(buf, "exact %s", eqns->vars.name[v]);
        const long len = strlen(eqns->vars.name[v]);
        if (eqns->vars.name[v][len + 1]) {
            strcat(buf, "-");
            strcat(buf, &eqns->vars.name[v][len + 1]);
        }
        name[v] = memory_strdup(buf);
    }
    equations_create_user(eqns, (void *)name, func, n_user);
    for (long v = 0; v < n_user; ++v) free(name[v]);
}

void equations_set_scalar(Equations *eqns, long scalar, double value)
{
    eqns->scalar.value[scalar] = value;
}

void equations_set_convective_flux(Equations *eqns, const char *name)
{
    strcpy(eqns->flux.name_conv, name);
    eqns->flux.conv = eqns->flux.select_conv(name);
}

void equations_set_viscous_flux(Equations *eqns, const char *name)
{
    strcpy(eqns->flux.name_visc, name);
    eqns->flux.visc = eqns->flux.select_visc(name);
}

void equations_set_prepare(Equations *eqns, Prepare *prepare)
{
    eqns->source.prepare = prepare;
}

void equations_set_source(Equations *eqns, Function *func)
{
    eqns->source.func = func;
}

void equations_set_limiter(Equations *eqns, const char *name, double k)
{
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const ALIAS(cv, eqns->mesh->cell.volume);
    memory_free(&eqns->limiter.eps2);
    if (!name) {
        strcpy(eqns->limiter.name, "");
        eqns->limiter.func = 0;
    }
    else if (!strcmp(name, "barth jespersen") || !strcmp(name, "minmod")) {
        strcpy(eqns->limiter.name, name);
        eqns->limiter.func = barth_jespersen;
        eqns->limiter.eps2 = memory_calloc(n_inner_cells, sizeof(*eqns->limiter.eps2));
    }
    else if (!strcmp(name, "venkatakrishnan") || !strcmp(name, "venk")) {
        strcpy(eqns->limiter.name, name);
        eqns->limiter.func = venkatakrishnan;
        eqns->limiter.eps2 = memory_calloc(n_inner_cells, sizeof(*eqns->limiter.eps2));

        eqns->limiter.k = k;
        for (long i = 0; i < n_inner_cells; ++i)
            eqns->limiter.eps2[i] = pow(k * pow(cv[i], 1.0 / N_DIMS), 3);
    }
    else
        error("unsupported limiter function '%s'", name);
}

void equations_set_initial_condition(Equations *eqns, Function *func)
{
    const long n_vars = eqns->n_vars;
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const ALIAS(x, eqns->mesh->cell.center);
    ALIAS(update, eqns->boundary);
    double(*u)[n_vars] = (void *)eqns->vars.u;

    for (long i = 0; i < n_inner_cells; ++i) {
        func(u[i], 0, x[i], 0);
        update(eqns, u[i]);
    }
}

void equations_set_initial_state(Equations *eqns, const double *state)
{
    const long n_vars = eqns->n_vars;
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    ALIAS(update, eqns->boundary);
    double(*u)[n_vars] = (void *)eqns->vars.u;

    for (long i = 0; i < n_inner_cells; ++i) {
        for (long v = 0; v < n_vars; ++v) u[i][v] = state[v];
        update(eqns, u[i]);
    }
}

void equations_set_boundary_condition(Equations *eqns, const char *entity, const char *bc,
                                      const double *state, Function *func)
{
    const long n_entities = eqns->mesh->n_entities;
    const ALIAS(name, eqns->mesh->entity.name);
    for (long e = 0; e < n_entities; ++e) {
        if (strcmp(name[e], entity)) continue;
        strcpy(eqns->bc.name[e], bc);
        eqns->bc.state[e] = state;
        if (strcmp(bc, "custom"))
            eqns->bc.apply[e] = eqns->bc.select(bc);
        else if (!func)
            error("custom boundary condition for entity '%s' requires function", entity);
        eqns->bc.func[e] = func;
        return;
    }
    error("could not find entity '%s'", entity);
}

static void compute_output_dimensions(char (*name)[NAMELEN], long **dim, long n_names)
{
    (*dim) = memory_calloc(n_names, sizeof(*(*dim)));
    for (long n = 0, i = 0; i < n_names; ++i) {
        (*dim)[n] += 1;
        char *dash = strrchr(name[i], '-');
        if (dash) *dash = 0;
        if (i == n_names - 1 || strncmp(name[i + 1], name[i], strlen(name[i]))) n += 1;
    }
}
