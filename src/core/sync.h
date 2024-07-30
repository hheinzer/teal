#pragma once

#include "equations.h"
#include "mesh.h"

/* Collectively compute the smallest value of 'v' among all ranks. */
#define sync_min(v) _Generic(v, long: x__sync_min_long, double: x__sync_min_double)(v)
long x__sync_min_long(long v);
double x__sync_min_double(double v);

/* Collectively compute the largest value of 'v' among all ranks. */
#define sync_max(v) _Generic(v, long: x__sync_max_long, double: x__sync_max_double)(v)
long x__sync_max_long(long v);
double x__sync_max_double(double v);

/* Collectively compute the sum of 'v' among all ranks. */
#define sync_sum(v) _Generic(v, long: x__sync_sum_long, double: x__sync_sum_double)(v)
long x__sync_sum_long(long v);
double x__sync_sum_double(double v);

/* Collectively compute the sum of the exclusive scan of 'v' on all ranks. */
#define sync_exscan_sum(v) \
    _Generic(v, long: x__sync_exscan_sum_long, double: x__sync_exscan_sum_double)(v)
long x__sync_exscan_sum_long(long v);
double x__sync_exscan_sum_double(double v);

/* Synchronize 'nu' variables of buffer 'u' among all ranks using the sync info of 'mesh'. */
void sync_all(const Mesh *mesh, double *u, long nu);

/* Begin the synchronization of 'nu' variables of buffer 'u' among all ranks using the sync info of
 * 'eqns'. This function must be followed by calls to sync_wait and sync_end. */
void sync_begin(Equations *eqns, double *u, long nu);

/* Wait until all receiving operations have completed using the sync info of 'eqns'. */
void sync_wait(Equations *eqns);

/* Wait until all sending operations have completed using the sync info of 'eqns'. */
void sync_end(Equations *eqns);
