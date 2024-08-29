#include "rotate.h"

#include "euler.h"

void euler_global_to_local(const Matrix3d b, const double *g, double *l)
{
    l[D] = g[D];
    l[U] = b[X][X] * g[U] + b[X][Y] * g[V] + b[X][Z] * g[W];
    l[V] = b[Y][X] * g[U] + b[Y][Y] * g[V] + b[Y][Z] * g[W];
    l[W] = b[Z][X] * g[U] + b[Z][Y] * g[V] + b[Z][Z] * g[W];
    l[P] = g[P];
}

void euler_local_to_global(const Matrix3d b, const double *l, double *g)
{
    g[D] = l[D];
    g[DU] = b[X][X] * l[DU] + b[Y][X] * l[DV] + b[Z][X] * l[DW];
    g[DV] = b[X][Y] * l[DU] + b[Y][Y] * l[DV] + b[Z][Y] * l[DW];
    g[DW] = b[X][Z] * l[DU] + b[Y][Z] * l[DV] + b[Z][Z] * l[DW];
    g[DE] = l[DE];
}
