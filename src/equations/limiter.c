#include "limiter.h"

#include "teal/isclose.h"
#include "teal/utils.h"

double barth_jespersen(const Vector3d *dx, const Vector3d dudx, double u, double umin, double umax,
                       double, long n)
{
    // Blazek 2015, sec. 5.3.5
    double psi = 1;
    for (long i = 0; i < n; ++i) {
        double delta2 = 0;
        for (long d = 0; d < N_DIMS; ++d) delta2 += dudx[d] * dx[i][d];
        if (is_close(delta2, 0)) continue;
        const double delta1 = (delta2 > 0 ? umax - u : umin - u);
        psi = min(psi, delta1 / delta2);
    }
    return psi;
}

double venkatakrishnan(const Vector3d *dx, const Vector3d dudx, double u, double umin, double umax,
                       double eps2, long n)
{
    // Blazek 2015, sec. 5.3.5
    double psi = 1;
    for (long i = 0; i < n; ++i) {
        double delta2 = 0;
        for (long d = 0; d < N_DIMS; ++d) delta2 += dudx[d] * dx[i][d];
        if (is_close(delta2, 0)) continue;
        const double delta1 = (delta2 > 0 ? umax - u : umin - u);
        psi = min(psi, ((sq(delta1) + eps2) * delta2 + 2 * sq(delta2) * delta1) /
                           (sq(delta1) + 2 * sq(delta2) + delta1 * delta2 + eps2) / delta2);
    }
    return psi;
}

void equations_limiter(Equations *eqns)
{
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const long n_vars = eqns->n_vars;
    const alias(i_cell, eqns->mesh->cell.i_cell);
    const alias(cell, eqns->mesh->cell.cell);
    const alias(dx, eqns->mesh->cell.to_cell);
    const alias(eps2, eqns->limiter.eps2);
    alias(limiter, eqns->limiter.limiter);
    const double(*u)[n_vars] = (void *)eqns->vars.u;
    Vector3d(*dudx)[n_vars] = (void *)eqns->vars.dudx;

    for (long i = 0; i < n_inner_cells; ++i) {
        for (long v = 0; v < n_vars; ++v) {
            double umin = u[i][v];
            double umax = u[i][v];
            for (long j = i_cell[i]; j < i_cell[i + 1]; ++j) {
                umin = min(umin, u[cell[j]][v]);
                umax = max(umax, u[cell[j]][v]);
            }
            const long n = i_cell[i + 1] - i_cell[i];
            const double psi = limiter(&dx[i_cell[i]], dudx[i][v], u[i][v], umin, umax, eps2[i], n);
            for (long d = 0; d < N_DIMS; ++d) dudx[i][v][d] *= psi;
        }
    }
}
