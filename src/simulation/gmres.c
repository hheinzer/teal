#include "gmres.h"

#include <assert.h>
#include <math.h>
#include <string.h>

#include "equations.h"
#include "sync.h"
#include "teal.h"
#include "utils.h"

static void matvec(const Equations *eqns, const void *conserved_, const void *derivative_,
                   const void *basis_, void *result_, double time, double step, double fd_scale)
{
    int num_cells = eqns->mesh->cells.num;
    int num_inner = eqns->mesh->cells.num_inner;
    int stride_p = eqns->primitive.stride;
    int stride_c = eqns->conserved.stride;
    assert(num_cells > 0 && num_inner > 0 && stride_p > 0 && stride_c > 0);
    double *property = eqns->properties.data;
    Convert *cons2prim = eqns->primitive.func.convert;

    const double (*conserved)[stride_c] = conserved_;
    const double (*derivative)[stride_c] = derivative_;
    const double (*basis)[stride_c] = basis_;
    double (*result)[stride_c] = result_;

    double eps = fd_scale / sync_norm(*basis, num_inner * stride_c);

    double perturbed[stride_c];
    double (*primitive1)[stride_p] = teal_calloc(num_cells, (int)sizeof(*primitive1));
    for (int i = 0; i < num_inner; i++) {
        for (int j = 0; j < stride_c; j++) {
            perturbed[j] = conserved[i][j] + (eps * basis[i][j]);
        }
        cons2prim(primitive1[i], perturbed, property);
    }

    double (*derivative1)[stride_c] = teal_calloc(num_inner, (int)sizeof(*derivative1));
    equations_derivative(eqns, primitive1, derivative1, time);

    for (int i = 0; i < num_inner; i++) {
        for (int j = 0; j < stride_c; j++) {
            result[i][j] = basis[i][j] - (step * (derivative1[i][j] - derivative[i][j]) / eps);
        }
    }

    teal_free(derivative1);
    teal_free(primitive1);
}

void gmres(const Equations *eqns, const void *conserved_, const void *derivative_,
           const void *residual_, void *increment_, double time, double step, double norm,
           double fd_scale, const NewtonKrylov *ctx)
{
    double tol = ctx->krylov_tolerance;
    int dim = ctx->krylov_dimension;

    int num_inner = eqns->mesh->cells.num_inner;
    int stride_c = eqns->conserved.stride;

    const double (*conserved)[stride_c] = conserved_;
    const double (*derivative)[stride_c] = derivative_;
    const double (*residual)[stride_c] = residual_;
    double (*increment)[stride_c] = increment_;

    memset(increment, 0, num_inner * sizeof(*increment));
    if (norm <= 0) {
        return;
    }

    double *rhs = teal_calloc(dim + 1, (int)sizeof(*rhs));
    rhs[0] = norm;

    double (*basis)[num_inner][stride_c] = teal_calloc(dim + 1, (int)sizeof(*basis));
    for (int i = 0; i < num_inner; i++) {
        for (int j = 0; j < stride_c; j++) {
            basis[0][i][j] = -residual[i][j] / norm;
        }
    }

    double tol_norm = tol * norm;

    double (*result)[stride_c] = teal_calloc(num_inner, (int)sizeof(*result));
    double (*hess)[dim] = teal_calloc(dim + 1, (int)sizeof(*hess));
    double *cosine = teal_calloc(dim, (int)sizeof(*cosine));
    double *sine = teal_calloc(dim, (int)sizeof(*sine));
    int iter = 0;
    for (; iter < dim; iter++) {
        matvec(eqns, conserved, derivative, basis[iter], result, time, step, fd_scale);

        for (int i = 0; i < iter + 1; i++) {
            hess[i][iter] = sync_dot(*basis[i], *result, num_inner * stride_c);
            for (int j = 0; j < num_inner; j++) {
                for (int k = 0; k < stride_c; k++) {
                    result[j][k] -= hess[i][iter] * basis[i][j][k];
                }
            }
        }
        hess[iter + 1][iter] = sync_norm(*result, num_inner * stride_c);

        for (int i = 0; i < iter; i++) {
            double tmp = (cosine[i] * hess[i][iter]) + (sine[i] * hess[i + 1][iter]);
            hess[i + 1][iter] = (-sine[i] * hess[i][iter]) + (cosine[i] * hess[i + 1][iter]);
            hess[i][iter] = tmp;
        }

        double rot = sqrt(sq(hess[iter][iter]) + sq(hess[iter + 1][iter]));
        cosine[iter] = hess[iter][iter] / rot;
        sine[iter] = hess[iter + 1][iter] / rot;
        hess[iter][iter] = rot;

        rhs[iter + 1] = -sine[iter] * rhs[iter];
        rhs[iter] = cosine[iter] * rhs[iter];

        if (fabs(rhs[iter + 1]) < tol_norm) {
            iter += 1;
            break;
        }

        for (int i = 0; i < num_inner; i++) {
            for (int j = 0; j < stride_c; j++) {
                basis[iter + 1][i][j] = result[i][j] / hess[iter + 1][iter];
            }
        }
    }

    double *coef = teal_calloc(dim, (int)sizeof(*coef));
    for (int i = iter - 1; i >= 0; i--) {
        coef[i] = rhs[i];
        for (int j = i + 1; j < iter; j++) {
            coef[i] -= hess[i][j] * coef[j];
        }
        coef[i] /= hess[i][i];
        for (int j = 0; j < num_inner; j++) {
            for (int k = 0; k < stride_c; k++) {
                increment[j][k] += coef[i] * basis[i][j][k];
            }
        }
    }

    teal_free(coef);
    teal_free(sine);
    teal_free(cosine);
    teal_free(hess);
    teal_free(result);
    teal_free(basis);
    teal_free(rhs);
}
