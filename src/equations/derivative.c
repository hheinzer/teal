#include <assert.h>
#include <math.h>
#include <string.h>

#include "equations.h"
#include "sync.h"
#include "teal/arena.h"
#include "teal/utils.h"

static void integrate_convective_flux(const Equations *eqns, void *variable_, void *derivative_)
{
    Arena save = arena_save();

    int num_faces = eqns->mesh->faces.num;
    int num_inner = eqns->mesh->faces.num_inner;
    int off_ghost = eqns->mesh->faces.off_ghost;
    Adjacent *cell = eqns->mesh->faces.cell;
    scalar *area = eqns->mesh->faces.area;
    matrix *basis = eqns->mesh->faces.basis;

    int len = eqns->variables.len;
    int stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Convective *convective = eqns->convective.flux;

    scalar(*variable)[stride] = variable_;
    scalar(*derivative)[len] = derivative_;
    scalar *flux = arena_malloc(len, sizeof(*flux));

    Request req = sync_variables(eqns, variable, stride);

    for (int i = 0; i < num_inner; i++) {
        int left = cell[i].left;
        int right = cell[i].right;
        convective(flux, variable[left], variable[right], property, &basis[i]);
        for (int j = 0; j < len; j++) {
            derivative[left][j] += flux[j] * area[i];
            derivative[right][j] -= flux[j] * area[i];
        }
    }
    for (int i = num_inner; i < off_ghost; i++) {
        int left = cell[i].left;
        int right = cell[i].right;
        convective(flux, variable[left], variable[right], property, &basis[i]);
        for (int j = 0; j < len; j++) {
            derivative[left][j] += flux[j] * area[i];
        }
    }

    sync_wait(eqns, req.recv);

    for (int i = off_ghost; i < num_faces; i++) {
        int left = cell[i].left;
        int right = cell[i].right;
        convective(flux, variable[left], variable[right], property, &basis[i]);
        for (int j = 0; j < len; j++) {
            derivative[left][j] += flux[j] * area[i];
        }
    }

    sync_wait(eqns, req.send);

    arena_load(save);
}

static void reconstruct(scalar *variable_k, const scalar *variable, const vector *gradient,
                        const vector *offset, int stride)
{
    for (int i = 0; i < stride; i++) {
        variable_k[i] = variable[i] + ((gradient[i].x * offset->x) + (gradient[i].y * offset->y) +
                                       (gradient[i].z * offset->z));
    }
}

static void integrate_reconstructed_convective_flux(const Equations *eqns, void *variable_,
                                                    void *derivative_, void *gradient_)
{
    Arena save = arena_save();

    int num_faces = eqns->mesh->faces.num;
    int num_inner = eqns->mesh->faces.num_inner;
    int off_ghost = eqns->mesh->faces.off_ghost;
    Adjacent *cell = eqns->mesh->faces.cell;
    scalar *area = eqns->mesh->faces.area;
    matrix *basis = eqns->mesh->faces.basis;
    Offset *offset = eqns->mesh->faces.offset;

    int len = eqns->variables.len;
    int stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Convective *convective = eqns->convective.flux;

    scalar(*variable)[stride] = variable_;
    scalar(*derivative)[len] = derivative_;
    vector(*gradient)[stride] = gradient_;
    scalar *variable_l = arena_malloc(stride, sizeof(*variable_l));
    scalar *variable_r = arena_malloc(stride, sizeof(*variable_r));
    scalar *flux = arena_malloc(len, sizeof(*flux));

    Request req = sync_gradients(eqns, gradient, stride);

    for (int i = 0; i < num_inner; i++) {
        int left = cell[i].left;
        int right = cell[i].right;
        reconstruct(variable_l, variable[left], gradient[left], &offset[i].left, stride);
        reconstruct(variable_r, variable[right], gradient[right], &offset[i].right, stride);
        convective(flux, variable_l, variable_r, property, &basis[i]);
        for (int j = 0; j < len; j++) {
            derivative[left][j] += flux[j] * area[i];
            derivative[right][j] -= flux[j] * area[i];
        }
    }
    for (int i = num_inner; i < off_ghost; i++) {
        int left = cell[i].left;
        int right = cell[i].right;
        reconstruct(variable_l, variable[left], gradient[left], &offset[i].left, stride);
        convective(flux, variable_l, variable[right], property, &basis[i]);
        for (int j = 0; j < len; j++) {
            derivative[left][j] += flux[j] * area[i];
        }
    }

    sync_wait(eqns, req.recv);

    for (int i = off_ghost; i < num_faces; i++) {
        int left = cell[i].left;
        int right = cell[i].right;
        reconstruct(variable_l, variable[left], gradient[left], &offset[i].left, stride);
        reconstruct(variable_r, variable[right], gradient[right], &offset[i].right, stride);
        convective(flux, variable_l, variable_r, property, &basis[i]);
        for (int j = 0; j < len; j++) {
            derivative[left][j] += flux[j] * area[i];
        }
    }

    sync_wait(eqns, req.send);

    arena_load(save);
}

static void evaluate_source(const Equations *eqns, void *variable_, void *derivative_, scalar time)
{
    Arena save = arena_save();

    int num = eqns->mesh->cells.num_inner;
    vector *center = eqns->mesh->cells.center;

    int len = eqns->variables.len;
    int stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Source *compute = eqns->source;

    scalar(*variable)[stride] = variable_;
    scalar(*derivative)[len] = derivative_;
    scalar *source = arena_malloc(len, sizeof(*source));

    for (int i = 0; i < num; i++) {
        compute(source, variable[i], property, center[i], time);
        for (int j = 0; j < len; j++) {
            derivative[i][j] += source[j];
        }
    }

    arena_load(save);
}

void *equations_derivative(const Equations *eqns, void *variable_, void *derivative_, scalar time)
{
    assert(eqns && variable_ && isfinite(time) && time >= 0);

    int num = eqns->mesh->cells.num_inner;
    int len = eqns->variables.len;

    scalar(*derivative)[len] = derivative_;
    if (!derivative) {
        derivative = arena_malloc(num, sizeof(*derivative));
    }
    memset(derivative, 0, num * sizeof(*derivative));

    equations_boundary(eqns, variable_, time);
    switch (eqns->space_order) {
        case 1: {
            integrate_convective_flux(eqns, variable_, derivative);
            break;
        }
        case 2: {
            void *gradient_ = equations_gradient(eqns, variable_);
            if (eqns->limiter.compute) {
                equations_limiter(eqns, variable_, gradient_);
            }
            if (eqns->convective.flux && eqns->viscous.flux) {
                error("TODO");
            }
            else if (eqns->convective.flux) {
                integrate_reconstructed_convective_flux(eqns, variable_, derivative, gradient_);
            }
            else if (eqns->viscous.flux) {
                error("TODO");
            }
            break;
        }
        default: error("invalid space order -- '%d'", eqns->space_order);
    }

    scalar *volume = eqns->mesh->cells.volume;
    for (int i = 0; i < num; i++) {
        for (int j = 0; j < len; j++) {
            derivative[i][j] = -derivative[i][j] / volume[i];
        }
    }

    if (eqns->source) {
        evaluate_source(eqns, variable_, derivative, time);
    }

    return derivative;
}
