#pragma once

#include <mpi.h>

extern struct Sync {
    MPI_Comm comm;
    int rank;
    int size;
} sync;

// Initialize MPI and cache rank and size.
void sync_init(int *argc, char ***argv);

// Finalize MPI and release the cached communicator.
void sync_deinit(void);

// Replace the cached communicator and refresh rank and size.
void sync_reinit(MPI_Comm comm);

long sync_lmin(long val);
long sync_lmax(long val);
long sync_lsum(long val);

double sync_fmin(double val);
double sync_fmax(double val);
double sync_fsum(double val);

// Return the exclusive prefix sum of a long value.
long sync_exsum(long val);

// Return the global dot product of two arrays.
double sync_dot(const double *lhs, const double *rhs, long num);

// Return the global Euclidean norm of an array.
double sync_norm(const double *arr, long num);
