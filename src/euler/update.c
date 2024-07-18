#include "update.h"

#include "core/utils.h"
#include "euler.h"

void prim_to_cons(const Equations *eqns, double *u)
{
    const double gamma = eqns->scalar.value[GAMMA];
    u[DU] = u[D] * u[U];
    u[DV] = u[D] * u[V];
    u[DW] = u[D] * u[W];
    u[DE] = u[P] / (gamma - 1) + 0.5 * (u[DU] * u[U] + u[DV] * u[V] + u[DW] * u[W]);
}

void cons_to_prim(const Equations *eqns, double *u)
{
    const double gamma = eqns->scalar.value[GAMMA];
    u[D] = max(EPS, u[D]);
    u[U] = u[DU] / u[D];
    u[V] = u[DV] / u[D];
    u[W] = u[DW] / u[D];
    u[P] = max(EPS, (gamma - 1) * (u[DE] - 0.5 * (u[DU] * u[U] + u[DV] * u[V] + u[DW] * u[W])));
}
