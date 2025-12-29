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

long sync_all_lmin(long val);
long sync_all_lmax(long val);
long sync_all_lsum(long val);

double sync_all_fmin(double val);
double sync_all_fmax(double val);
double sync_all_fsum(double val);

long sync_exsum(long val);
