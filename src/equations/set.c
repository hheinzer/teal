#include <assert.h>
#include <string.h>

#include "equations.h"
#include "limiter.h"
#include "teal.h"

void equations_set_space_order(Equations *eqns, int space_order)
{
    assert(eqns && 1 <= space_order && space_order <= 2);
    eqns->space_order = space_order;
}

void equations_set_convective_flux(Equations *eqns, const char *name)
{
    assert(eqns);

    if (!name) {
        eqns->convective.name[0] = 0;
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
        eqns->viscous.name[0] = 0;
        eqns->viscous.flux = 0;
        return;
    }

    assert(eqns->viscous.select);
    strcpy(eqns->viscous.name, name);
    eqns->viscous.flux = eqns->viscous.select(name);
}

void equations_set_limiter(Equations *eqns, const char *name, double coefficient)
{
    assert(eqns);

    teal_free(eqns->limiter.parameter);
    eqns->limiter.parameter = 0;

    if (!name) {
        eqns->limiter.name[0] = 0;
        eqns->limiter.kind = NONE;
        return;
    }

    strcpy(eqns->limiter.name, name);

    if (!strcmp(name, "barth jespersen") || !strcmp(name, "minmod")) {
        eqns->limiter.kind = BARTH_JESPERSEN;
    }
    else if (!strcmp(name, "venkatakrishnan")) {
        eqns->limiter.kind = VENKATAKRISHNAN;
        eqns->limiter.parameter = venkatakrishnan_parameter(eqns, coefficient);
    }
    else {
        teal_error("invalid limiter (%s)", name);
    }
}

void equations_set_boundary_condition(Equations *eqns, const char *entity, const char *name,
                                      Compute *compute, const void *context)
{
    assert(eqns && entity);

    int num_inner = eqns->mesh->entities.num_inner;

    for (int i = 0; i < eqns->boundary.num; i++) {
        if (!strcmp(eqns->mesh->entities.name[num_inner + i], entity)) {
            if (!name || !strcmp(name, "custom")) {
                assert(compute);
                strcpy(eqns->boundary.name[i], "custom");
                eqns->boundary.compute[i] = compute;
                eqns->boundary.condition[i] = 0;
            }
            else {
                assert(eqns->boundary.select);
                strcpy(eqns->boundary.name[i], name);
                eqns->boundary.compute[i] = 0;
                eqns->boundary.condition[i] = eqns->boundary.select(name);
            }
            eqns->boundary.context[i] = context;
            return;
        }
    }
    teal_error("invalid entity (%s)", entity);
}

void equations_set_initial_condition(Equations *eqns, const char *entity, Compute *compute,
                                     const void *initial)
{
    assert(eqns && entity);

    Vector *center = eqns->mesh->cells.center;

    int num_inner = eqns->mesh->entities.num_inner;
    int *cell_off = eqns->mesh->entities.cell_off;

    int stride = eqns->primitive.stride;
    double (*primitive)[stride] = eqns->primitive.data;
    double *property = eqns->properties.data;

    for (int i = 0; i < num_inner; i++) {
        if (!strcmp(eqns->mesh->entities.name[i], entity)) {
            for (int j = cell_off[i]; j < cell_off[i + 1]; j++) {
                if (compute) {
                    assert(!initial);
                    compute(primitive[j], property, center[j], 0, 0);
                }
                else {
                    assert(initial);
                    memcpy(primitive[j], initial, sizeof(*primitive));
                }
            }
            return;
        }
    }
    teal_error("invalid entity (%s)", entity);
}

void equations_set_property(Equations *eqns, const char *name, double property)
{
    assert(eqns && name);
    for (int i = 0; i < eqns->properties.num; i++) {
        if (!strcmp(eqns->properties.name[i], name)) {
            eqns->properties.data[i] = property;
            return;
        }
    }
    teal_error("invalid property (%s)", name);
}

void equations_set_source(Equations *eqns, Compute *compute)
{
    assert(eqns && compute);
    eqns->source.compute = compute;
}
