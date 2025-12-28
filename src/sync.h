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
