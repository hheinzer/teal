#include "rotate.h"

#include "euler.h"

void global_to_local(const double *n, const double *g, double *l)
{
    const double(*r)[N_DIMS] = (const void *)n;
    l[D] = g[D];
    l[U] = r[X][X] * g[U] + r[X][Y] * g[V] + r[X][Z] * g[W];
    l[V] = r[Y][X] * g[U] + r[Y][Y] * g[V] + r[Y][Z] * g[W];
    l[W] = r[Z][X] * g[U] + r[Z][Y] * g[V] + r[Z][Z] * g[W];
    l[P] = g[P];
}

void local_to_global(const double *n, const double *l, double *g)
{
    const double(*r)[N_DIMS] = (const void *)n;
    g[D] = l[D];
    g[DU] = r[X][X] * l[DU] + r[Y][X] * l[DV] + r[Z][X] * l[DW];
    g[DV] = r[X][Y] * l[DU] + r[Y][Y] * l[DV] + r[Z][Y] * l[DW];
    g[DW] = r[X][Z] * l[DU] + r[Y][Z] * l[DV] + r[Z][Z] * l[DW];
    g[DE] = l[DE];
}
