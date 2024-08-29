#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "equations.h"
#include "limiter.h"
#include "mesh/find.h"
#include "teal/utils.h"

void equations_set_timestep(Equations *eqns, TimeStep *compute)
{
    eqns->timestep.compute = compute;
}

void equations_set_update(Equations *eqns, Update *boundary, Update *advance)
{
    eqns->update.boundary = boundary;
    eqns->update.advance = advance;
}

void equations_set_select_convective_flux(Equations *eqns, SelectConvFlux *select)
{
    eqns->conv.select = select;
}

void equations_set_select_viscous_flux(Equations *eqns, SelectViscFlux *select)
{
    eqns->visc.select = select;
}

void equations_set_select_boundary_condition(Equations *eqns, SelectApplyBC *select)
{
    eqns->bc.select = select;
}

void equations_set_scalar(Equations *eqns, long name, double value)
{
    eqns->scalar.value[name] = value;
}

void equations_set_convective_flux(Equations *eqns, const char *name)
{
    strcpy(eqns->conv.name, name);
    eqns->conv.flux = eqns->conv.select(name);
}

void equations_set_viscous_flux(Equations *eqns, const char *name)
{
    strcpy(eqns->visc.name, name);
    eqns->visc.flux = eqns->visc.select(name);
}

void equations_set_source(Equations *eqns, Compute *compute, Prepare *prepare)
{
    eqns->source.compute = compute;
    eqns->source.prepare = prepare;
}

void equations_set_space_order(Equations *eqns, long space_order)
{
    assert(1 <= space_order && space_order <= 2);
    eqns->space_order = space_order;
}

void equations_set_limiter(Equations *eqns, const char *name, double k)
{
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const alias(cv, eqns->mesh->cell.volume);
    alias(eps2, eqns->limiter.eps2);

    strcpy(eqns->limiter.name, name);
    if (!strcmp(name, "none"))
        eqns->limiter.limiter = 0;
    else if (!strcmp(name, "barth jespersen") || !strcmp(name, "minmod"))
        eqns->limiter.limiter = barth_jespersen;
    else if (!strcmp(name, "venkatakrishnan") || !strcmp(name, "venk")) {
        eqns->limiter.limiter = venkatakrishnan;
        eqns->limiter.k = k;
        for (long i = 0; i < n_inner_cells; ++i) eps2[i] = pow(k * pow(cv[i], 1.0 / N_DIMS), 3);
    }
    else
        abort();
}

void equations_set_boundary_condition(Equations *eqns, const char *entity, const char *name,
                                      const double *state, Compute *compute)
{
    const long E = mesh_find_entity(eqns->mesh, entity);
    strcpy(eqns->bc.name[E], name);
    eqns->bc.state[E] = state;
    if (strcmp(name, "custom"))
        eqns->bc.apply[E] = eqns->bc.select(name);
    else
        assert(compute);
    eqns->bc.compute[E] = compute;
}

void equations_set_initial_condition(Equations *eqns, Compute *compute)
{
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const long n_vars = eqns->n_vars;
    const alias(x, eqns->mesh->cell.center);
    alias(update, eqns->update.boundary);
    double(*u)[n_vars] = (void *)eqns->vars.u;

    for (long i = 0; i < n_inner_cells; ++i) {
        compute(u[i], 0, x[i], 0);
        update(eqns, u[i]);
    }
}

void equations_set_initial_state(Equations *eqns, const double *state)
{
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const long n_vars = eqns->n_vars;
    alias(update, eqns->update.boundary);
    double(*u)[n_vars] = (void *)eqns->vars.u;

    for (long i = 0; i < n_inner_cells; ++i) {
        for (long v = 0; v < n_vars; ++v) u[i][v] = state[v];
        update(eqns, u[i]);
    }
}
