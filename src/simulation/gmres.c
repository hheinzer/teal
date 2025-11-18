#include "gmres.h"

#include <math.h>
#include <string.h>

#include "teal/arena.h"
#include "teal/sync.h"
#include "teal/utils.h"

static void matvec(const Equations *eqns, const void *variable_, const void *derivative_,
                   const void *basis_, void *result_, scalar time, scalar step, scalar fd_scale)
{
    Arena save = arena_save();

    int num_cells = eqns->mesh->cells.num;
    int num = eqns->mesh->cells.num_inner;
    int len = eqns->variables.len;
    int stride = eqns->variables.stride;
    scalar *property = eqns->properties.data;
    Update *primitive = eqns->variables.primitive;

    const scalar(*variable)[stride] = variable_;
    const scalar(*derivative)[len] = derivative_;
    const scalar(*basis)[len] = basis_;
    scalar(*result)[len] = result_;

    scalar eps = fd_scale / sync_fnorm(*basis, num * len);

    scalar(*variable1)[stride] = arena_malloc(num_cells, sizeof(*variable1));
    for (int i = 0; i < num; i++) {
        for (int j = 0; j < len; j++) {
            variable1[i][j] = variable[i][j] + (eps * basis[i][j]);
        }
        primitive(variable1[i], property);
    }

    scalar(*derivative1)[len] = equations_derivative(eqns, variable1, 0, time);

    for (int i = 0; i < num; i++) {
        for (int j = 0; j < len; j++) {
            result[i][j] = basis[i][j] - (step * (derivative1[i][j] - derivative[i][j]) / eps);
        }
    }

    arena_load(save);
}

void gmres(const Equations *eqns, const void *variable_, const void *derivative_,
           const void *residual_, void *increment_, scalar time, scalar step, scalar norm,
           scalar fd_scale, const NewtonKrylov *ctx)
{
    Arena save = arena_save();

    scalar tol = ctx->krylov_tolerance;
    int dim = ctx->krylov_dimension;

    int num = eqns->mesh->cells.num_inner;
    int len = eqns->variables.len;
    int stride = eqns->variables.stride;

    const scalar(*variable)[stride] = variable_;
    const scalar(*derivative)[len] = derivative_;
    const scalar(*residual)[len] = residual_;
    scalar(*increment)[len] = increment_;

    memset(increment, 0, num * sizeof(*increment));
    if (norm <= 0) {
        return;
    }

    scalar *rhs = arena_malloc(dim + 1, sizeof(*rhs));
    rhs[0] = norm;

    scalar(*basis)[num][len] = arena_malloc(dim + 1, sizeof(*basis));
    for (int i = 0; i < num; i++) {
        for (int j = 0; j < len; j++) {
            basis[0][i][j] = -residual[i][j] / norm;
        }
    }

    scalar tol_norm = tol * norm;

    scalar(*result)[len] = arena_malloc(num, sizeof(*result));
    scalar(*hess)[dim] = arena_malloc(dim + 1, sizeof(*hess));
    scalar *cosine = arena_malloc(dim, sizeof(*cosine));
    scalar *sine = arena_malloc(dim, sizeof(*sine));
    int iter = 0;
    for (; iter < dim; iter++) {
        matvec(eqns, variable, derivative, basis[iter], result, time, step, fd_scale);

        for (int i = 0; i < iter + 1; i++) {
            hess[i][iter] = sync_fdot(*basis[i], *result, num * len);
            for (int j = 0; j < num; j++) {
                for (int k = 0; k < len; k++) {
                    result[j][k] -= hess[i][iter] * basis[i][j][k];
                }
            }
        }
        hess[iter + 1][iter] = sync_fnorm(*result, num * len);

        for (int i = 0; i < iter; i++) {
            scalar tmp = (cosine[i] * hess[i][iter]) + (sine[i] * hess[i + 1][iter]);
            hess[i + 1][iter] = (-sine[i] * hess[i][iter]) + (cosine[i] * hess[i + 1][iter]);
            hess[i][iter] = tmp;
        }

        scalar rot = sqrt(pow2(hess[iter][iter]) + pow2(hess[iter + 1][iter]));
        cosine[iter] = hess[iter][iter] / rot;
        sine[iter] = hess[iter + 1][iter] / rot;
        hess[iter][iter] = rot;

        rhs[iter + 1] = -sine[iter] * rhs[iter];
        rhs[iter] = cosine[iter] * rhs[iter];

        if (fabs(rhs[iter + 1]) < tol_norm) {
            iter += 1;
            break;
        }

        for (int i = 0; i < num; i++) {
            for (int j = 0; j < len; j++) {
                basis[iter + 1][i][j] = result[i][j] / hess[iter + 1][iter];
            }
        }
    }

    scalar *coeff = arena_malloc(dim, sizeof(*coeff));
    for (int i = iter - 1; i >= 0; i--) {
        coeff[i] = rhs[i];
        for (int j = i + 1; j < iter; j++) {
            coeff[i] -= hess[i][j] * coeff[j];
        }
        coeff[i] /= hess[i][i];
        for (int j = 0; j < num; j++) {
            for (int k = 0; k < len; k++) {
                increment[j][k] += coeff[i] * basis[i][j][k];
            }
        }
    }

    arena_load(save);
}
