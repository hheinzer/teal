#pragma once

#include <mpi.h>

extern struct Sync {
    MPI_Comm comm;
    int rank, size;
} sync;

void sync_init(int *argc, char ***argv);

#define sync_min(v) _Generic(v, long: x__sync_min_long, double: x__sync_min_double)(v)
long x__sync_min_long(long v);
double x__sync_min_double(double v);

#define sync_max(v) _Generic(v, long: x__sync_max_long, double: x__sync_max_double)(v)
long x__sync_max_long(long v);
double x__sync_max_double(double v);

#define sync_sum(v) _Generic(v, long: x__sync_sum_long, double: x__sync_sum_double)(v)
long x__sync_sum_long(long v);
double x__sync_sum_double(double v);

#define sync_exsum(v) _Generic(v, long: x__sync_exsum_long, double: x__sync_exsum_double)(v)
long x__sync_exsum_long(long v);
double x__sync_exsum_double(double v);

double sync_dot(const double *a, const double *b, long n);

double sync_norm(const double *a, long n);

#define sync_gatherv(a, b, n, root) \
    _Generic(*a, long: x__sync_gatherv_long, double: x__sync_gatherv_double)(a, b, n, root)
void x__sync_gatherv_long(long *a, const long *b, long n, int root);
void x__sync_gatherv_double(double *a, const double *b, long n, int root);

#define sync_scatterv(a, b, n, root) \
    _Generic(*a, long: x__sync_scatterv_long, double: x__sync_scatterv_double)(a, b, n, root)
void x__sync_scatterv_long(long *a, const long *b, long n, int root);
void x__sync_scatterv_double(double *a, const double *b, long n, int root);

void sync_irecv(const long *j_recv, MPI_Request *req, double *u, long ldu);

void sync_isend(const long *i_send, const long *send, const double *u, MPI_Request *req,
                double *buf, long ldu);

void sync_waitall(MPI_Request *req);

void sync_finalize(void);
