#ifndef SYNC_H
#define SYNC_H

#include "equations.h"

long sync_min_long(const long v);
long *sync_min_long_n(long *v, const long n);
double sync_min_double(const double v);
double *sync_min_double_n(double *v, const long n);
#define sync_min(v, ...)         \
    _Generic(v,                  \
        long: sync_min_long,     \
        long *: sync_min_long_n, \
        double: sync_min_double, \
        double *: sync_min_double_n)(v __VA_OPT__(, __VA_ARGS__))

long sync_max_long(const long v);
long *sync_max_long_n(long *v, const long n);
double sync_max_double(const double v);
double *sync_max_double_n(double *v, const long n);
#define sync_max(v, ...)         \
    _Generic(v,                  \
        long: sync_max_long,     \
        long *: sync_max_long_n, \
        double: sync_max_double, \
        double *: sync_max_double_n)(v __VA_OPT__(, __VA_ARGS__))

long sync_sum_long(const long v);
long *sync_sum_long_n(long *v, const long n);
double sync_sum_double(const double v);
double *sync_sum_double_n(double *v, const long n);
#define sync_sum(v, ...)         \
    _Generic(v,                  \
        long: sync_sum_long,     \
        long *: sync_sum_long_n, \
        double: sync_sum_double, \
        double *: sync_sum_double_n)(v __VA_OPT__(, __VA_ARGS__))

long sync_exscan_sum_long(const long v);
long *sync_exscan_sum_long_n(long *v, const long n);
double sync_exscan_sum_double(const double v);
double *sync_exscan_sum_double_n(double *v, const long n);
#define sync_exscan_sum(v, ...)         \
    _Generic(v,                         \
        long: sync_exscan_sum_long,     \
        long *: sync_exscan_sum_long_n, \
        double: sync_exscan_sum_double, \
        double *: sync_exscan_sum_double_n)(v __VA_OPT__(, __VA_ARGS__))

void sync_all(const Mesh *mesh, const long nu, double (*u)[nu]);

void sync_u_begin(Equations *eqns);

void sync_u_wait(Equations *eqns);

void sync_u_end(Equations *eqns);

void sync_dudx_begin(Equations *eqns);

void sync_dudx_wait(Equations *eqns);

void sync_dudx_end(Equations *eqns);

#endif
