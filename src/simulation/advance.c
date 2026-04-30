#include "advance.h"

#include <math.h>

#include "equations.h"
#include "gmres.h"
#include "sync.h"
#include "teal.h"
#include "utils.h"

double euler(const Equations *eqns, double *time, void *residual, double max_step, double courant,
             const void *context)
{
    (void)context;

    int num_inner = eqns->mesh->cells.num_inner;
    int stride_p = eqns->primitive.stride;
    int stride_c = eqns->conserved.stride;
    double *property = eqns->properties.data;
    Convert *prim2cons = eqns->conserved.func.convert;
    Convert *cons2prim = eqns->primitive.func.convert;

    double (*primitive)[stride_p] = eqns->primitive.data;
    double (*conserved)[stride_c] = eqns->conserved.data;

    double step0 = courant * equations_timestep(eqns, primitive);
    double step = fmin(step0, max_step);

    double (*derivative)[stride_c] = teal_calloc(num_inner, (int)sizeof(*derivative));
    equations_derivative(eqns, primitive, derivative, *time);

    for (int i = 0; i < num_inner; i++) {
        prim2cons(conserved[i], primitive[i], property);
        for (int j = 0; j < stride_c; j++) {
            conserved[i][j] += derivative[i][j] * step;
        }
        cons2prim(primitive[i], conserved[i], property);
    }

    *time += step;
    equations_residual(eqns, derivative, residual);

    teal_free(derivative);
    return step0;
}

double lserk(const Equations *eqns, double *time, void *residual, double max_step, double courant,
             const void *context_)
{
    static const double alpha_[3][6][7] = {
        {
            {0.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.3333, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.1481, 0.4000, 1.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.0833, 0.2069, 0.4265, 1.0000, 0.0000, 0.0000},
            {0.0000, 0.0533, 0.1263, 0.2375, 0.4414, 1.0000, 0.0000},
            {0.0000, 0.0370, 0.0851, 0.1521, 0.2562, 0.4512, 1.0000},
        },
        {
            {0.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.4242, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.1918, 0.4929, 1.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.1084, 0.2602, 0.5052, 1.0000, 0.0000, 0.0000},
            {0.0000, 0.0695, 0.1602, 0.2898, 0.5060, 1.0000, 0.0000},
            {0.0000, 0.0482, 0.1085, 0.1885, 0.3050, 0.5063, 1.0000},
        },
        {
            {0.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.6612, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.2884, 0.5010, 1.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.1666, 0.3027, 0.5275, 1.0000, 0.0000, 0.0000},
            {0.0000, 0.1067, 0.1979, 0.3232, 0.5201, 1.0000, 0.0000},
            {0.0000, 0.0742, 0.1393, 0.2198, 0.3302, 0.5181, 1.0000},
        },
    };

    const RungeKutta *context = context_;
    int time_order = context->time_order;
    int num_stages = context->num_stages;
    const double *alpha = alpha_[time_order - 1][num_stages - 1];

    int num_inner = eqns->mesh->cells.num_inner;
    int stride_p = eqns->primitive.stride;
    int stride_c = eqns->conserved.stride;
    double *property = eqns->properties.data;
    Convert *prim2cons = eqns->conserved.func.convert;
    Convert *cons2prim = eqns->primitive.func.convert;

    double (*primitive)[stride_p] = eqns->primitive.data;
    double (*conserved)[stride_c] = eqns->conserved.data;

    double step0 = courant * equations_timestep(eqns, primitive);
    double step = fmin(step0, max_step);

    double (*conserved0)[stride_c] = teal_calloc(num_inner, (int)sizeof(*conserved0));
    for (int i = 0; i < num_inner; i++) {
        prim2cons(conserved0[i], primitive[i], property);
    }

    double (*derivative)[stride_c] = teal_calloc(num_inner, (int)sizeof(*derivative));
    for (int i = 0; i < num_stages; i++) {
        equations_derivative(eqns, primitive, derivative, *time + (alpha[i] * step));
        for (int j = 0; j < num_inner; j++) {
            for (int k = 0; k < stride_c; k++) {
                conserved[j][k] = conserved0[j][k] + (derivative[j][k] * alpha[i + 1] * step);
            }
            cons2prim(primitive[j], conserved[j], property);
        }
    }

    *time += step;
    equations_residual(eqns, derivative, residual);

    teal_free(conserved0);
    teal_free(derivative);
    return step0;
}

double implicit_euler(const Equations *eqns, double *time, void *residual, double max_step,
                      double courant, const void *context_)
{
    const NewtonKrylov *context = context_;
    double tol = context->newton_tolerance;

    int num_inner = eqns->mesh->cells.num_inner;
    int stride_p = eqns->primitive.stride;
    int stride_c = eqns->conserved.stride;
    double *property = eqns->properties.data;
    Convert *prim2cons = eqns->conserved.func.convert;
    Convert *cons2prim = eqns->primitive.func.convert;

    double (*primitive)[stride_p] = eqns->primitive.data;
    double (*conserved)[stride_c] = eqns->conserved.data;

    double step0 = courant * equations_timestep(eqns, primitive);
    double step = fmin(step0, max_step);

    double (*conserved0)[stride_c] = teal_calloc(num_inner, (int)sizeof(*conserved0));
    for (int i = 0; i < num_inner; i++) {
        prim2cons(conserved[i], primitive[i], property);
        copy(conserved0[i], conserved[i], stride_c, (int)sizeof(double));
    }

    *time += step;
    double (*derivative)[stride_c] = teal_calloc(num_inner, (int)sizeof(*derivative));
    equations_derivative(eqns, primitive, derivative, *time);

    double (*newton_residual)[stride_c] = teal_calloc(num_inner, (int)sizeof(*newton_residual));
    for (int i = 0; i < num_inner; i++) {
        for (int j = 0; j < stride_c; j++) {
            newton_residual[i][j] = -derivative[i][j] * step;
        }
    }

    double norm = sync_norm(*newton_residual, num_inner * stride_c);

    static const double fd_relative = 1e-6;
    double cons_norm = sync_norm(*conserved, num_inner * stride_c);
    double fd_scale = fd_relative * fmax(1, cons_norm) / fmax(1, norm);

    double tol_norm = tol * norm;
    int max_iter = 128 * (num_inner + 1);

    double (*increment)[stride_c] = teal_calloc(num_inner, (int)sizeof(*increment));
    for (int iter = 0; iter < max_iter && norm > tol_norm; iter++) {
        gmres(eqns, conserved, derivative, newton_residual, increment, *time, step, norm, fd_scale,
              context);
        for (int i = 0; i < num_inner; i++) {
            for (int j = 0; j < stride_c; j++) {
                conserved[i][j] += increment[i][j];
            }
            cons2prim(primitive[i], conserved[i], property);
        }

        equations_derivative(eqns, primitive, derivative, *time);
        for (int i = 0; i < num_inner; i++) {
            for (int j = 0; j < stride_c; j++) {
                newton_residual[i][j] =
                    conserved[i][j] - conserved0[i][j] - (derivative[i][j] * step);
            }
        }

        norm = sync_norm(*newton_residual, num_inner * stride_c);
        cons_norm = sync_norm(*conserved, num_inner * stride_c);
        fd_scale = fd_relative * fmax(1, cons_norm) / fmax(1, norm);
    }

    equations_residual(eqns, derivative, residual);

    teal_free(increment);
    teal_free(newton_residual);
    teal_free(derivative);
    teal_free(conserved0);
    return step0;
}
