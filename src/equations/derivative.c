#include <math.h>

#include "equations.h"
#include "sync.h"
#include "teal/arena.h"
#include "teal/assert.h"
#include "teal/utils.h"

static void integrate_convective_flux(const Equations *eqns, void *variable_, void *derivative_)
{
    Arena save = arena_save();

    number num_faces = eqns->mesh->faces.num;
    number num_inner = eqns->mesh->faces.num_inner;
    number off_ghost = eqns->mesh->faces.off_ghost;
    Adjacent *cell = eqns->mesh->faces.cell;
    scalar *area = eqns->mesh->faces.area;
    matrix *basis = eqns->mesh->faces.basis;

    number len = eqns->variables.len;
    number stride = eqns->variables.stride;
    scalar *property = eqns->properties.property;
    Convective *convective = eqns->convective.flux;

    scalar(*variable)[stride] = variable_;
    scalar(*derivative)[len] = derivative_;
    scalar *flux = arena_malloc(len, sizeof(*flux));

    Request req = sync_variables(eqns, variable, stride);

    for (number i = 0; i < num_inner; i++) {
        number left = cell[i].left;
        number right = cell[i].right;
        convective(flux, variable[left], variable[right], property, basis[i]);
        for (number j = 0; j < len; j++) {
            derivative[left][j] += flux[j] * area[i];
            derivative[right][j] -= flux[j] * area[i];
        }
    }
    for (number i = num_inner; i < off_ghost; i++) {
        number left = cell[i].left;
        number right = cell[i].right;
        convective(flux, variable[left], variable[right], property, basis[i]);
        for (number j = 0; j < len; j++) {
            derivative[left][j] += flux[j] * area[i];
        }
    }

    sync_wait(eqns, req.recv);

    for (number i = off_ghost; i < num_faces; i++) {
        number left = cell[i].left;
        number right = cell[i].right;
        convective(flux, variable[left], variable[right], property, basis[i]);
        for (number j = 0; j < len; j++) {
            derivative[left][j] += flux[j] * area[i];
        }
    }

    sync_wait(eqns, req.send);

    arena_load(save);
}

static void evaluate_source(const Equations *eqns, void *variable_, void *derivative_, scalar time)
{
    Arena save = arena_save();

    number num = eqns->mesh->cells.num_inner;
    vector *center = eqns->mesh->cells.center;

    number len = eqns->variables.len;
    number stride = eqns->variables.stride;
    scalar *property = eqns->properties.property;
    Source *compute = eqns->source;

    scalar(*variable)[stride] = variable_;
    scalar(*derivative)[len] = derivative_;
    scalar *source = arena_malloc(len, sizeof(*source));

    for (number i = 0; i < num; i++) {
        compute(source, variable[i], property, center[i], time);
        for (number j = 0; j < len; j++) {
            derivative[i][j] += source[j];
        }
    }

    arena_load(save);
}

void *equations_derivative(const Equations *eqns, void *variable_, scalar time)
{
    assert(eqns && variable_ && isfinite(time) && time >= 0);

    number num = eqns->mesh->cells.num_inner;
    number len = eqns->variables.len;
    scalar(*derivative)[len] = arena_calloc(num, sizeof(*derivative));

    switch (eqns->space_order) {
        case 1: integrate_convective_flux(eqns, variable_, derivative); break;
        default: error("invalid space order -- '%td'", eqns->space_order);
    }

    scalar *volume = eqns->mesh->cells.volume;
    for (number i = 0; i < num; i++) {
        for (number j = 0; j < len; j++) {
            derivative[i][j] = -derivative[i][j] / volume[i];
        }
    }

    if (eqns->source) {
        evaluate_source(eqns, variable_, derivative, time);
    }

    return derivative;
}
