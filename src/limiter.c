#include "limiter.h"

#include "global.h"
#include "mesh.h"
#include "utils.h"

double limiter_barth_jespersen(const double, const double u, const double umin, const double umax,
                               const double *dudx, const double (*dx)[N_DIMS], const long nj)
{
    // Blazek 2015, sec. 5.3.5
    double psi = 1;
    for (long j = 0; j < nj; ++j) {
        double delta2 = EPS;
        for (long d = 0; d < N_DIMS; ++d) delta2 += dudx[d] * dx[j][d];
        const double du = (delta2 > 0 ? umax - u : umin - u);
        psi = MIN(psi, du / delta2);
    }
    return psi;
}

double limiter_venkatakrishnan(const double eps2, const double u, const double umin,
                               const double umax, const double *dudx, const double (*dx)[N_DIMS],
                               const long nj)
{
    // Blazek 2015, sec. 5.3.5
    double psi = 1;
    for (long j = 0; j < nj; ++j) {
        double delta2 = EPS;
        for (long d = 0; d < N_DIMS; ++d) delta2 += dudx[d] * dx[j][d];
        const double du = (delta2 > 0 ? umax - u : umin - u);
        psi = MIN(psi, ((du * du + eps2) * delta2 + 2 * delta2 * delta2 * du) /
                           (du * du + 2 * delta2 * delta2 + du * delta2 + eps2) / delta2);
    }
    return psi;
}
