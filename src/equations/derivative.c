#include <assert.h>
#include <string.h>

#include "equations.h"
#include "exchange.h"
#include "limiter.h"
#include "teal.h"

static void initialize_derivative(const Equations *eqns, void *derivative_)
{
    int num_inner = eqns->mesh->cells.num_inner;
    int stride = eqns->conserved.stride;

    double (*derivative)[stride] = derivative_;
    memset(derivative, 0, num_inner * sizeof(*derivative));
}

static void integrate_convective_flux_O1(const Equations *eqns, void *primitive_, void *derivative_)
{
    int num_faces = eqns->mesh->faces.num;
    int num_inner = eqns->mesh->faces.num_inner;
    int off_boundary = eqns->mesh->faces.off_boundary;
    Pair *cell_idx = eqns->mesh->faces.cell_idx;
    double *area = eqns->mesh->faces.area;
    Matrix *basis = eqns->mesh->faces.basis;

    int stride_p = eqns->primitive.stride;
    int stride_c = eqns->conserved.stride;
    double (*primitive)[stride_p] = primitive_;
    double (*derivative)[stride_c] = derivative_;
    double *property = eqns->properties.data;

    Convective *convective = eqns->convective.flux;
    double flux[stride_c];

    Exchange exchange = equations_exchange(eqns, primitive, stride_p);

    for (int i = 0; i < num_inner; i++) {
        int left = cell_idx[i].left;
        int right = cell_idx[i].right;
        convective(flux, primitive[left], primitive[right], property, basis[i]);
        for (int j = 0; j < stride_c; j++) {
            derivative[left][j] += flux[j] * area[i];
            derivative[right][j] -= flux[j] * area[i];
        }
    }

    for (int i = num_inner; i < off_boundary; i++) {
        int left = cell_idx[i].left;
        int right = cell_idx[i].right;
        convective(flux, primitive[left], primitive[right], property, basis[i]);
        for (int j = 0; j < stride_c; j++) {
            derivative[left][j] += flux[j] * area[i];
        }
    }

    equations_exchange_wait_recv(eqns, exchange);

    for (int i = off_boundary; i < num_faces; i++) {
        int left = cell_idx[i].left;
        int right = cell_idx[i].right;
        convective(flux, primitive[left], primitive[right], property, basis[i]);
        for (int j = 0; j < stride_c; j++) {
            derivative[left][j] += flux[j] * area[i];
        }
    }

    equations_exchange_wait_send(eqns, exchange);
}

static void average(void *avg_, const void *lhs_, const void *rhs_, int len)
{
    double *avg = avg_;
    const double *lhs = lhs_;
    const double *rhs = rhs_;
    for (int i = 0; i < len; i++) {
        avg[i] = (lhs[i] + rhs[i]) / 2;
    }
}

static void correct(Vector *gradient_m, const double *primitive_l, const double *primitive_r,
                    VectorNorm correction, int stride)
{
    for (int i = 0; i < stride; i++) {
        double projection = vector_dot(gradient_m[i], correction.unit);
        double adjust = projection - ((primitive_r[i] - primitive_l[i]) / correction.norm);
        vector_isub(&gradient_m[i], vector_mul(adjust, correction.unit));
    }
}

static void integrate_viscous_flux_O2(const Equations *eqns, const void *primitive_,
                                      void *derivative_, void *gradient_)
{
    int num_faces = eqns->mesh->faces.num;
    int num_inner = eqns->mesh->faces.num_inner;
    int off_boundary = eqns->mesh->faces.off_boundary;
    Pair *cell_idx = eqns->mesh->faces.cell_idx;
    double *area = eqns->mesh->faces.area;
    Matrix *basis = eqns->mesh->faces.basis;
    VectorNorm *correction = eqns->mesh->faces.correction;

    int stride_p = eqns->primitive.stride;
    int stride_c = eqns->conserved.stride;
    const double (*primitive)[stride_p] = primitive_;
    double (*derivative)[stride_c] = derivative_;
    Vector(*gradient)[stride_p] = gradient_;
    double *property = eqns->properties.data;

    Viscous *viscous = eqns->viscous.flux;
    double primitive_m[stride_p];
    Vector gradient_m[stride_p];
    double flux[stride_c];

    Exchange exchange = equations_exchange(eqns, gradient, 3 * stride_p);

    for (int i = 0; i < num_inner; i++) {
        int left = cell_idx[i].left;
        int right = cell_idx[i].right;
        average(primitive_m, primitive[left], primitive[right], stride_p);
        average(gradient_m, gradient[left], gradient[right], 3 * stride_p);
        correct(gradient_m, primitive[left], primitive[right], correction[i], stride_p);
        viscous(flux, primitive_m, gradient_m, property, basis[i]);
        for (int j = 0; j < stride_c; j++) {
            derivative[left][j] -= flux[j] * area[i];
            derivative[right][j] += flux[j] * area[i];
        }
    }

    for (int i = num_inner; i < off_boundary; i++) {
        int left = cell_idx[i].left;
        int right = cell_idx[i].right;
        average(primitive_m, primitive[left], primitive[right], stride_p);
        memcpy(gradient_m, gradient[left], sizeof(gradient_m));
        correct(gradient_m, primitive[left], primitive[right], correction[i], stride_p);
        viscous(flux, primitive_m, gradient_m, property, basis[i]);
        for (int j = 0; j < stride_c; j++) {
            derivative[left][j] -= flux[j] * area[i];
        }
    }

    equations_exchange_wait_recv(eqns, exchange);

    for (int i = off_boundary; i < num_faces; i++) {
        int left = cell_idx[i].left;
        int right = cell_idx[i].right;
        average(primitive_m, primitive[left], primitive[right], stride_p);
        average(gradient_m, gradient[left], gradient[right], 3 * stride_p);
        correct(gradient_m, primitive[left], primitive[right], correction[i], stride_p);
        viscous(flux, primitive_m, gradient_m, property, basis[i]);
        for (int j = 0; j < stride_c; j++) {
            derivative[left][j] -= flux[j] * area[i];
        }
    }

    equations_exchange_wait_send(eqns, exchange);
}

static void reconstruct(double *primitive_k, const double *primitive, const Vector *gradient,
                        Vector offset, int stride)
{
    for (int i = 0; i < stride; i++) {
        primitive_k[i] = primitive[i] + vector_dot(gradient[i], offset);
    }
}

static void integrate_convective_flux_O2(const Equations *eqns, const void *primitive_,
                                         void *derivative_, void *gradient_)
{
    int num_faces = eqns->mesh->faces.num;
    int num_inner = eqns->mesh->faces.num_inner;
    int off_boundary = eqns->mesh->faces.off_boundary;
    Pair *cell_idx = eqns->mesh->faces.cell_idx;
    double *area = eqns->mesh->faces.area;
    Matrix *basis = eqns->mesh->faces.basis;
    VectorPair *offset = eqns->mesh->faces.offset;

    int stride_p = eqns->primitive.stride;
    int stride_c = eqns->conserved.stride;
    const double (*primitive)[stride_p] = primitive_;
    double (*derivative)[stride_c] = derivative_;
    Vector(*gradient)[stride_p] = gradient_;
    double *property = eqns->properties.data;

    Convective *convective = eqns->convective.flux;
    double primitive_l[stride_p];
    double primitive_r[stride_p];
    double flux[stride_c];

    Exchange exchange = equations_exchange(eqns, gradient, 3 * stride_p);

    for (int i = 0; i < num_inner; i++) {
        int left = cell_idx[i].left;
        int right = cell_idx[i].right;
        reconstruct(primitive_l, primitive[left], gradient[left], offset[i].left, stride_p);
        reconstruct(primitive_r, primitive[right], gradient[right], offset[i].right, stride_p);
        convective(flux, primitive_l, primitive_r, property, basis[i]);
        for (int j = 0; j < stride_c; j++) {
            derivative[left][j] += flux[j] * area[i];
            derivative[right][j] -= flux[j] * area[i];
        }
    }

    for (int i = num_inner; i < off_boundary; i++) {
        int left = cell_idx[i].left;
        int right = cell_idx[i].right;
        reconstruct(primitive_l, primitive[left], gradient[left], offset[i].left, stride_p);
        convective(flux, primitive_l, primitive[right], property, basis[i]);
        for (int j = 0; j < stride_c; j++) {
            derivative[left][j] += flux[j] * area[i];
        }
    }

    equations_exchange_wait_recv(eqns, exchange);

    for (int i = off_boundary; i < num_faces; i++) {
        int left = cell_idx[i].left;
        int right = cell_idx[i].right;
        reconstruct(primitive_l, primitive[left], gradient[left], offset[i].left, stride_p);
        reconstruct(primitive_r, primitive[right], gradient[right], offset[i].right, stride_p);
        convective(flux, primitive_l, primitive_r, property, basis[i]);
        for (int j = 0; j < stride_c; j++) {
            derivative[left][j] += flux[j] * area[i];
        }
    }

    equations_exchange_wait_send(eqns, exchange);
}

static void finalize_derivative(const Equations *eqns, const void *primitive_, void *derivative_,
                                double time)
{
    int num_inner = eqns->mesh->cells.num_inner;
    double *volume = eqns->mesh->cells.volume;
    Vector *center = eqns->mesh->cells.center;

    int stride_p = eqns->primitive.stride;
    int stride_c = eqns->conserved.stride;
    const double (*primitive)[stride_p] = primitive_;
    double (*derivative)[stride_c] = derivative_;
    double *property = eqns->properties.data;

    Compute *compute = eqns->source.compute;
    double source[stride_c];

    if (!compute) {
        for (int i = 0; i < num_inner; i++) {
            for (int j = 0; j < stride_c; j++) {
                derivative[i][j] = -derivative[i][j] / volume[i];
            }
        }
        return;
    }

    for (int i = 0; i < num_inner; i++) {
        compute(source, property, center[i], time, primitive[i]);
        for (int j = 0; j < stride_c; j++) {
            derivative[i][j] = (-derivative[i][j] / volume[i]) + source[j];
        }
    }
}

void equations_derivative(const Equations *eqns, void *primitive, void *derivative, double time)
{
    assert(eqns && primitive && derivative);

    initialize_derivative(eqns, derivative);

    equations_boundary(eqns, primitive, time);

    switch (eqns->space_order) {
        case 1: integrate_convective_flux_O1(eqns, primitive, derivative); break;
        case 2: {
            int num_cells = eqns->mesh->cells.num;
            int stride = eqns->primitive.stride;

            Vector(*gradient)[stride] = teal_calloc(num_cells, sizeof(*gradient));
            equations_gradient(eqns, primitive, gradient);

            if (eqns->viscous.flux) {
                integrate_viscous_flux_O2(eqns, primitive, derivative, gradient);
            }

            if (eqns->limiter.kind != NONE) {
                equations_limiter(eqns, primitive, gradient);
            }

            if (eqns->convective.flux) {
                integrate_convective_flux_O2(eqns, primitive, derivative, gradient);
            }

            teal_free(gradient);
            break;
        }
        default: teal_error("invalid space order (%d)", eqns->space_order);
    }

    finalize_derivative(eqns, primitive, derivative, time);
}
