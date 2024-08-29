#include "update.h"

#include "euler.h"
#include "teal/utils.h"

void euler_prim_to_cons(const Equations *eqns, double *u)
{
    const double gamma = eqns->scalar.value[GAMMA];
    u[DU] = u[D] * u[U];
    u[DV] = u[D] * u[V];
    u[DW] = u[D] * u[W];
    u[DE] = u[P] / (gamma - 1) + 0.5 * (u[DU] * u[U] + u[DV] * u[V] + u[DW] * u[W]);
}

void euler_cons_to_prim(const Equations *eqns, double *u)
{
    const double gamma = eqns->scalar.value[GAMMA];
    u[D] = max(1e-8, u[D]);
    u[U] = u[DU] / u[D];
    u[V] = u[DV] / u[D];
    u[W] = u[DW] / u[D];
    u[P] = max(1e-8, (gamma - 1) * (u[DE] - 0.5 * (u[DU] * u[U] + u[DV] * u[V] + u[DW] * u[W])));
}
