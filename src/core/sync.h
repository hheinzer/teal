#ifndef SYNC_H
#define SYNC_H

#include "equations.h"
#include "mesh.h"

long sync_min_long(const long v);
double sync_min_double(const double v);
#define sync_min(v) _Generic(v, long: sync_min_long, double: sync_min_double)(v)

long sync_max_long(const long v);
double sync_max_double(const double v);
#define sync_max(v) _Generic(v, long: sync_max_long, double: sync_max_double)(v)

long sync_sum_long(const long v);
double sync_sum_double(const double v);
#define sync_sum(v) _Generic(v, long: sync_sum_long, double: sync_sum_double)(v)

long sync_exscan_sum_long(const long v);
double sync_exscan_sum_double(const double v);
#define sync_exscan_sum(v) \
    _Generic(v, long: sync_exscan_sum_long, double: sync_exscan_sum_double)(v)

long *sync_gather_long(long *v, const long n);
double *sync_gather_double(double *v, const long n);
#define sync_gather(v, n) _Generic(*v, long: sync_gather_long, double: sync_gather_double)(v, n)

void sync_all(const Mesh *mesh, const long nu, double (*u)[nu]);

void sync_begin(Equations *eqns, const long nu, double (*u)[nu]);

void sync_wait(Equations *eqns);

void sync_end(Equations *eqns);

#endif
