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

// Return the global minimum of an int value.
int sync_min(int val);

// Return the global maximum of an int value.
int sync_max(int val);

// Return the global sum of an int value.
int sync_sum(int val);

// Return the global minimum of a double value.
double sync_fmin(double val);

// Return the global maximum of a double value.
double sync_fmax(double val);

// Return the global sum of a double value.
double sync_fsum(double val);

// Return the exclusive prefix sum of an int value.
int sync_exsum(int val);

// Return the global dot product of two arrays.
double sync_dot(const double *lhs, const double *rhs, int num);

// Return the global Euclidean norm of an array.
double sync_norm(const double *arr, int num);
