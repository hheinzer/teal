#ifndef SYNC_H
#define SYNC_H

#include "equations.h"
#include "mesh.h"

#define sync_min(v) _Generic(v, long: x__sync_min_long, double: x__sync_min_double)(v)
long x__sync_min_long(long v);
double x__sync_min_double(double v);

#define sync_max(v) _Generic(v, long: x__sync_max_long, double: x__sync_max_double)(v)
long x__sync_max_long(long v);
double x__sync_max_double(double v);

#define sync_sum(v) _Generic(v, long: x__sync_sum_long, double: x__sync_sum_double)(v)
long x__sync_sum_long(long v);
double x__sync_sum_double(double v);

#define sync_exscan_sum(v) \
    _Generic(v, long: x__sync_exscan_sum_long, double: x__sync_exscan_sum_double)(v)
long x__sync_exscan_sum_long(long v);
double x__sync_exscan_sum_double(double v);

void sync_all(const Mesh *mesh, double *u, long nu);

void sync_begin(Equations *eqns, double *u, long nu);

void sync_wait(Equations *eqns);

void sync_end(Equations *eqns);

#endif
