#include <assert.h>
#include <string.h>

#include "equations.h"
#include "limiter.h"
#include "teal/utils.h"

void equations_set_space_order(Equations *eqns, long space_order)
{
    assert(eqns);
    if (!(1 <= space_order && space_order <= 2)) {
        error("invalid space order (%ld)", space_order);
    }
    eqns->space_order = space_order;
}

void equations_set_time_step(Equations *eqns, TimeStep *time_step)
{
    assert(eqns && time_step);
    eqns->time_step = time_step;
}

void equations_set_boundary_select(Equations *eqns, BoundarySelect *select)
{
    assert(eqns && select);
    eqns->boundary.select = select;
}

void equations_set_convective_select(Equations *eqns, ConvectiveSelect *select)
{
    assert(eqns && select);
    eqns->convective.select = select;
}

void equations_set_viscous_select(Equations *eqns, ViscousSelect *select)
{
    assert(eqns && select);
    eqns->viscous.select = select;
}

void equations_set_limiter(Equations *eqns, const char *name, scalar parameter)
{
    assert(eqns);
    if (!name) {
        eqns->limiter.compute = 0;
        return;
    }
    strcpy(eqns->limiter.name, name);
    if (!strcmp(name, "vanleer")) {
        eqns->limiter.compute = vanleer;
        return;
    }
    if (!strcmp(name, "vanalbada1")) {
        eqns->limiter.compute = vanalbada1;
        return;
    }
    if (!strcmp(name, "vanalbada2")) {
        eqns->limiter.compute = vanalbada2;
        return;
    }
    if (!strcmp(name, "mc")) {
        eqns->limiter.compute = mc;
        return;
    }
    if (!strcmp(name, "koren")) {
        eqns->limiter.compute = koren;
        return;
    }
    if (!strcmp(name, "minmod")) {
        eqns->limiter.compute = minmod;
        return;
    }
    if (!strcmp(name, "superbee")) {
        eqns->limiter.compute = superbee;
        return;
    }
    if (!strcmp(name, "venkatakrishnan")) {
        eqns->limiter.parameter = venkatakrishnan_parameter(eqns, parameter);
        eqns->limiter.compute = venkatakrishnan;
        return;
    }
    error("invalid limiter (%s)", name);
}

void equations_set_boundary_condition(Equations *eqns, const char *entity, const char *name,
                                      const void *reference, Compute *custom)
{
    assert(eqns && entity && name);
    long num_inner = eqns->mesh->entities.num_inner;
    for (long i = 0; i < eqns->boundary.num; i++) {
        if (!strcmp(eqns->mesh->entities.name[num_inner + i], entity)) {
            strcpy(eqns->boundary.name[i], name);
            eqns->boundary.reference[i] = reference;
            eqns->boundary.custom[i] = custom;
            if (!strcmp(name, "custom") && !custom) {
                error("missing custom function (%s)", entity);
            }
            else {
                assert(eqns->boundary.select);
                eqns->boundary.condition[i] = eqns->boundary.select(name);
            }
            return;
        }
    }
    error("invalid entity (%s)", entity);
}

void equations_set_convective_flux(Equations *eqns, const char *name)
{
    assert(eqns);
    if (!name) {
        eqns->convective.flux = 0;
        return;
    }
    assert(eqns->convective.select);
    strcpy(eqns->convective.name, name);
    eqns->convective.flux = eqns->convective.select(name);
}

void equations_set_viscous_flux(Equations *eqns, const char *name)
{
    assert(eqns);
    if (!name) {
        eqns->viscous.flux = 0;
        return;
    }
    assert(eqns->viscous.select);
    strcpy(eqns->viscous.name, name);
    eqns->viscous.flux = eqns->viscous.select(name);
}

void equations_set_initial_condition(Equations *eqns, const char *entity, Compute *compute,
                                     scalar time)
{
    assert(eqns && entity && compute && time >= 0);

    vector *center = eqns->mesh->cells.center;

    long num_inner = eqns->mesh->entities.num_inner;
    Name *name = eqns->mesh->entities.name;
    long *cell_off = eqns->mesh->entities.cell_off;

    long stride = eqns->variables.stride;
    scalar(*variable)[stride] = eqns->variables.data;
    Update *conserved = eqns->variables.conserved;
    scalar *property = eqns->properties.data;

    for (long i = 0; i < num_inner; i++) {
        if (!strcmp(name[i], entity)) {
            for (long j = cell_off[i]; j < cell_off[i + 1]; j++) {
                compute(variable[j], property, center[j], time, 0);
                conserved(variable[j], property);
            }
            return;
        }
    }
    error("invalid entity (%s)", entity);
}

void equations_set_initial_state(Equations *eqns, const char *entity, const void *state)
{
    assert(eqns && entity && state);

    long num_inner = eqns->mesh->entities.num_inner;
    Name *name = eqns->mesh->entities.name;
    long *cell_off = eqns->mesh->entities.cell_off;

    long stride = eqns->variables.stride;
    scalar(*variable)[stride] = eqns->variables.data;
    Update *conserved = eqns->variables.conserved;
    scalar *property = eqns->properties.data;

    for (long i = 0; i < num_inner; i++) {
        if (!strcmp(name[i], entity)) {
            for (long j = cell_off[i]; j < cell_off[i + 1]; j++) {
                memcpy(variable[j], state, sizeof(*variable));
                conserved(variable[j], property);
            }
            return;
        }
    }
    error("invalid entity (%s)", entity);
}

void equations_set_property(Equations *eqns, const char *name, scalar property)
{
    assert(eqns && name);
    for (long i = 0; i < eqns->properties.num; i++) {
        if (!strcmp(eqns->properties.name[i], name)) {
            eqns->properties.data[i] = property;
            return;
        }
    }
    error("invalid property (%s)", name);
}

void equations_set_user_source(Equations *eqns, Source *source)
{
    assert(eqns && source);
    eqns->source = source;
}
