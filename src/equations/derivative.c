#include <assert.h>
#include <math.h>
#include <string.h>

#include "equations.h"
#include "sync.h"
#include "teal/arena.h"
#include "teal/utils.h"
#include "teal/vector.h"

// Allocate and zero the inner-cell derivative buffer.
static void initialize_derivative(const Equations *eqns, void **derivative_)
{
    long num_inner = eqns->mesh->cells.num_inner;
    long len = eqns->variables.len;
    scalar(*derivative)[len] = *derivative_;
    if (!derivative) {
        derivative = arena_malloc(num_inner, sizeof(*derivative));
    }
    *derivative_ = memset(derivative, 0, num_inner * sizeof(*derivative));
}

// Integrate first-order convective fluxes over all faces.
static void integrate_convective_flux_O1(const Equations *eqns, void *variable_, void *derivative_)
{
    Arena save = arena_save();

    long num_faces = eqns->mesh->faces.num;
    long num_inner = eqns->mesh->faces.num_inner;
    long off_ghost = eqns->mesh->faces.off_ghost;
    Adjacent *cell = eqns->mesh->faces.cell;
    scalar *area = eqns->mesh->faces.area;
    Basis *basis = eqns->mesh->faces.basis;

    long len = eqns->variables.len;
    long stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Convective *convective = eqns->convective.flux;

    scalar(*variable)[stride] = variable_;
    scalar(*derivative)[len] = derivative_;
    scalar *flux = arena_calloc(len, sizeof(*flux));

    Request req = sync_variables(eqns, variable, stride);

    for (long i = 0; i < num_inner; i++) {
        long left = cell[i].left;
        long right = cell[i].right;
        convective(flux, variable[left], variable[right], property, &basis[i]);
        for (long j = 0; j < len; j++) {
            derivative[left][j] += flux[j] * area[i];
            derivative[right][j] -= flux[j] * area[i];
        }
    }
    for (long i = num_inner; i < off_ghost; i++) {
        long left = cell[i].left;
        long right = cell[i].right;
        convective(flux, variable[left], variable[right], property, &basis[i]);
        for (long j = 0; j < len; j++) {
            derivative[left][j] += flux[j] * area[i];
        }
    }

    sync_wait(eqns, req.recv);

    for (long i = off_ghost; i < num_faces; i++) {
        long left = cell[i].left;
        long right = cell[i].right;
        convective(flux, variable[left], variable[right], property, &basis[i]);
        for (long j = 0; j < len; j++) {
            derivative[left][j] += flux[j] * area[i];
        }
    }

    sync_wait(eqns, req.send);

    arena_load(save);
}

// Average left/right states.
static void average_variable(scalar *variable_m, const scalar *variable_l, const scalar *variable_r,
                             long stride)
{
    for (long j = 0; j < stride; j++) {
        variable_m[j] = (variable_l[j] + variable_r[j]) / 2;
    }
}

// Average gradients between left and right cells.
static void average_gradient(vector *gradient_m, const vector *gradient_l, const vector *gradient_r,
                             long stride)
{
    for (long j = 0; j < stride; j++) {
        gradient_m[j].x = (gradient_l[j].x + gradient_r[j].x) / 2;
        gradient_m[j].y = (gradient_l[j].y + gradient_r[j].y) / 2;
        gradient_m[j].z = (gradient_l[j].z + gradient_r[j].z) / 2;
    }
}

// Adjust averaged gradients so their projection along the face offset matches the state jump.
static void correct_gradient(vector *gradient_m, const scalar *variable_l, const scalar *variable_r,
                             const Correction *correction, long stride)
{
    for (long j = 0; j < stride; j++) {
        scalar projection = vector_dot(gradient_m[j], correction->unit);
        scalar adjust = (projection - ((variable_r[j] - variable_l[j]) / correction->norm));
        gradient_m[j].x -= adjust * correction->unit.x;
        gradient_m[j].y -= adjust * correction->unit.y;
        gradient_m[j].z -= adjust * correction->unit.z;
    }
}

// Integrate second-order viscous fluxes with averaged state/gradients.
static void integrate_viscous_flux_O2(const Equations *eqns, void *variable_, void *derivative_,
                                      void *gradient_)
{
    Arena save = arena_save();

    long num_faces = eqns->mesh->faces.num;
    long num_inner = eqns->mesh->faces.num_inner;
    long off_ghost = eqns->mesh->faces.off_ghost;
    Adjacent *cell = eqns->mesh->faces.cell;
    scalar *area = eqns->mesh->faces.area;
    Basis *basis = eqns->mesh->faces.basis;
    Correction *correction = eqns->mesh->faces.correction;

    long len = eqns->variables.len;
    long stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Viscous *viscous = eqns->viscous.flux;

    scalar(*variable)[stride] = variable_;
    scalar(*derivative)[len] = derivative_;
    vector(*gradient)[stride] = gradient_;
    scalar *variable_m = arena_calloc(stride, sizeof(*variable_m));
    vector *gradient_m = arena_calloc(stride, sizeof(*gradient_m));
    scalar *flux = arena_calloc(len, sizeof(*flux));

    Request req = sync_gradients(eqns, gradient, stride);

    for (long i = 0; i < num_inner; i++) {
        long left = cell[i].left;
        long right = cell[i].right;
        average_variable(variable_m, variable[left], variable[right], stride);
        average_gradient(gradient_m, gradient[left], gradient[right], stride);
        correct_gradient(gradient_m, variable[left], variable[right], &correction[i], stride);
        viscous(flux, variable_m, gradient_m, property, &basis[i]);
        for (long j = 0; j < len; j++) {
            derivative[left][j] += -flux[j] * area[i];
            derivative[right][j] -= -flux[j] * area[i];
        }
    }
    for (long i = num_inner; i < off_ghost; i++) {
        long left = cell[i].left;
        long right = cell[i].right;
        average_variable(variable_m, variable[left], variable[right], stride);
        memcpy(gradient_m, gradient[left], stride * sizeof(*gradient_m));
        correct_gradient(gradient_m, variable[left], variable[right], &correction[i], stride);
        viscous(flux, variable_m, gradient_m, property, &basis[i]);
        for (long j = 0; j < len; j++) {
            derivative[left][j] += -flux[j] * area[i];
        }
    }

    sync_wait(eqns, req.recv);

    for (long i = off_ghost; i < num_faces; i++) {
        long left = cell[i].left;
        long right = cell[i].right;
        average_variable(variable_m, variable[left], variable[right], stride);
        average_gradient(gradient_m, gradient[left], gradient[right], stride);
        correct_gradient(gradient_m, variable[left], variable[right], &correction[i], stride);
        viscous(flux, variable_m, gradient_m, property, &basis[i]);
        for (long j = 0; j < len; j++) {
            derivative[left][j] += -flux[j] * area[i];
        }
    }

    sync_wait(eqns, req.send);

    arena_load(save);
}

// Linear reconstruction of a face state using gradients and offsets.
static void reconstruct(scalar *variable_k, const scalar *variable, const vector *gradient,
                        vector offset, long stride)
{
    for (long i = 0; i < stride; i++) {
        variable_k[i] = variable[i] + vector_dot(gradient[i], offset);
    }
}

// Integrate second-order convective fluxes with reconstructed states.
static void integrate_convective_flux_O2(const Equations *eqns, void *variable_, void *derivative_,
                                         void *gradient_)
{
    Arena save = arena_save();

    long num_faces = eqns->mesh->faces.num;
    long num_inner = eqns->mesh->faces.num_inner;
    long off_ghost = eqns->mesh->faces.off_ghost;
    Adjacent *cell = eqns->mesh->faces.cell;
    scalar *area = eqns->mesh->faces.area;
    Basis *basis = eqns->mesh->faces.basis;
    Offset *offset = eqns->mesh->faces.offset;

    long len = eqns->variables.len;
    long stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Convective *convective = eqns->convective.flux;

    scalar(*variable)[stride] = variable_;
    scalar(*derivative)[len] = derivative_;
    vector(*gradient)[stride] = gradient_;
    scalar *variable_l = arena_calloc(stride, sizeof(*variable_l));
    scalar *variable_r = arena_calloc(stride, sizeof(*variable_r));
    scalar *flux = arena_calloc(len, sizeof(*flux));

    Request req = sync_gradients(eqns, gradient, stride);

    for (long i = 0; i < num_inner; i++) {
        long left = cell[i].left;
        long right = cell[i].right;
        reconstruct(variable_l, variable[left], gradient[left], offset[i].left, stride);
        reconstruct(variable_r, variable[right], gradient[right], offset[i].right, stride);
        convective(flux, variable_l, variable_r, property, &basis[i]);
        for (long j = 0; j < len; j++) {
            derivative[left][j] += flux[j] * area[i];
            derivative[right][j] -= flux[j] * area[i];
        }
    }
    for (long i = num_inner; i < off_ghost; i++) {
        long left = cell[i].left;
        long right = cell[i].right;
        reconstruct(variable_l, variable[left], gradient[left], offset[i].left, stride);
        convective(flux, variable_l, variable[right], property, &basis[i]);
        for (long j = 0; j < len; j++) {
            derivative[left][j] += flux[j] * area[i];
        }
    }

    sync_wait(eqns, req.recv);

    for (long i = off_ghost; i < num_faces; i++) {
        long left = cell[i].left;
        long right = cell[i].right;
        reconstruct(variable_l, variable[left], gradient[left], offset[i].left, stride);
        reconstruct(variable_r, variable[right], gradient[right], offset[i].right, stride);
        convective(flux, variable_l, variable_r, property, &basis[i]);
        for (long j = 0; j < len; j++) {
            derivative[left][j] += flux[j] * area[i];
        }
    }

    sync_wait(eqns, req.send);

    arena_load(save);
}

// Divide flux sums by cell volume and add optional source terms.
static void finalize_derivative(const Equations *eqns, const void *variable_, void *derivative_,
                                scalar time)
{
    Arena save = arena_save();

    long num_inner = eqns->mesh->cells.num_inner;
    scalar *volume = eqns->mesh->cells.volume;

    long len = eqns->variables.len;
    scalar(*derivative)[len] = derivative_;

    if (eqns->source) {
        vector *center = eqns->mesh->cells.center;

        long stride = eqns->variables.stride;
        scalar *property = eqns->properties.data;
        Source *compute = eqns->source;

        const scalar(*variable)[stride] = variable_;
        scalar *source = arena_calloc(len, sizeof(*source));

        for (long i = 0; i < num_inner; i++) {
            compute(source, variable[i], property, center[i], time);
            for (long j = 0; j < len; j++) {
                derivative[i][j] = (-derivative[i][j] / volume[i]) + source[j];
            }
        }
    }
    else {
        for (long i = 0; i < num_inner; i++) {
            for (long j = 0; j < len; j++) {
                derivative[i][j] = (-derivative[i][j] / volume[i]);
            }
        }
    }

    arena_load(save);
}

void *equations_derivative(const Equations *eqns, void *variable_, void *derivative_, scalar time)
{
    assert(eqns && variable_ && isfinite(time) && time >= 0);
    initialize_derivative(eqns, &derivative_);
    equations_boundary(eqns, variable_, time);
    switch (eqns->space_order) {
        case 1: {
            integrate_convective_flux_O1(eqns, variable_, derivative_);
            break;
        }
        case 2: {
            void *gradient_ = equations_gradient(eqns, variable_);
            if (eqns->viscous.flux) {
                integrate_viscous_flux_O2(eqns, variable_, derivative_, gradient_);
            }
            if (eqns->limiter.compute) {
                equations_limiter(eqns, variable_, gradient_);
            }
            if (eqns->convective.flux) {
                integrate_convective_flux_O2(eqns, variable_, derivative_, gradient_);
            }
            break;
        }
        default: error("invalid space order (%ld)", eqns->space_order);
    }
    finalize_derivative(eqns, variable_, derivative_, time);
    return derivative_;
}
