#include <string.h>

#include "equations.h"
#include "limiter.h"
#include "teal/assert.h"
#include "teal/utils.h"

void equations_set_space_order(Equations *eqns, number space_order)
{
    assert(eqns && 1 <= space_order && space_order <= 2);
    eqns->space_order = space_order;
}

void equations_set_timestep(Equations *eqns, Timestep *timestep)
{
    assert(eqns && timestep);
    eqns->timestep = timestep;
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
    if (!strcmp(name, "minmod")) {
        eqns->limiter.compute = minmod;
        return;
    }
    if (!strcmp(name, "venkatakrishnan")) {
        eqns->limiter.parameter = venkatakrishnan_parameter(eqns, parameter);
        eqns->limiter.compute = venkatakrishnan;
        return;
    }
    error("invalid limiter -- '%s'", name);
}

void equations_set_boundary_condition(Equations *eqns, const char *entity, const char *name,
                                      const void *reference, Compute *custom)
{
    assert(eqns && entity && name);
    for (number i = 0; i < eqns->boundary.num; i++) {
        if (!strcmp(eqns->boundary.entity[i], entity)) {
            strcpy(eqns->boundary.name[i], name);
            eqns->boundary.reference[i] = reference;
            eqns->boundary.custom[i] = custom;
            if (!strcmp(name, "custom") && !custom) {
                error("missing custom function for entity -- '%s'", entity);
            }
            else {
                assert(eqns->boundary.select);
                eqns->boundary.condition[i] = eqns->boundary.select(name);
            }
            return;
        }
    }
    error("invalid entity -- '%s'", entity);
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

    number num = eqns->mesh->entities.num_inner;
    Name *name = eqns->mesh->entities.name;
    number *cell_off = eqns->mesh->entities.cell_off;

    number stride = eqns->variables.stride;
    scalar(*variable)[stride] = eqns->variables.data;
    Update *conserved = eqns->variables.conserved;
    scalar *property = eqns->properties.data;

    for (number i = 0; i < num; i++) {
        if (!strcmp(name[i], entity)) {
            for (number j = cell_off[i]; j < cell_off[i + 1]; j++) {
                compute(variable[j], property, center[j], time, 0);
                if (conserved) {
                    conserved(variable[j], property);
                }
            }
            return;
        }
    }
    error("invalid entity -- '%s'", entity);
}

void equations_set_initial_state(Equations *eqns, const char *entity, const void *state)
{
    assert(eqns && entity && state);

    number num = eqns->mesh->entities.num_inner;
    Name *name = eqns->mesh->entities.name;
    number *cell_off = eqns->mesh->entities.cell_off;

    number stride = eqns->variables.stride;
    scalar(*variable)[stride] = eqns->variables.data;
    Update *conserved = eqns->variables.conserved;
    scalar *property = eqns->properties.data;

    for (number i = 0; i < num; i++) {
        if (!strcmp(name[i], entity)) {
            for (number j = cell_off[i]; j < cell_off[i + 1]; j++) {
                memcpy(variable[j], state, sizeof(*variable));
                if (conserved) {
                    conserved(variable[j], property);
                }
            }
            return;
        }
    }
    error("invalid entity -- '%s'", entity);
}

void equations_set_property(Equations *eqns, const char *name, scalar property)
{
    assert(eqns && name);
    for (number i = 0; i < eqns->properties.num; i++) {
        if (!strcmp(eqns->properties.name[i], name)) {
            eqns->properties.data[i] = property;
            return;
        }
    }
    error("invalid property -- '%s'", name);
}

void equations_set_user_source(Equations *eqns, Source *source)
{
    assert(eqns && source);
    eqns->source = source;
}
