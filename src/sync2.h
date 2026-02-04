#pragma once

#include <mpi.h>

extern struct Sync2 {
    MPI_Comm comm;
    int rank;
    int size;
} sync2;

// Initialize MPI.
void sync2_init(int *argc, char ***argv);

// Finalize MPI.
void sync2_deinit(void);

// Replace the communicator and refresh rank/size.
void sync2_reinit(MPI_Comm comm);

// In-place global minimum.
void sync2_min(void *buf, int num, MPI_Datatype type);

// In-place global maximum.
void sync2_max(void *buf, int num, MPI_Datatype type);

// In-place global sum.
void sync2_sum(void *buf, int num, MPI_Datatype type);

// In-place exclusive prefix sum; rank 0 receives zeros.
void sync2_prefix(void *buf, int num, MPI_Datatype type);

// Resize a derived datatype to `extent` bytes.
MPI_Datatype sync2_resized(MPI_Datatype type, MPI_Aint extent);
