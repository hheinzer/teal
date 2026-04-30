#pragma once

#include <mpi.h>

extern struct sync {
    MPI_Comm comm;
    int rank;
    int size;
} sync;

// Initialize MPI.
void sync_init(int *argc, char ***argv);

// Finalize MPI.
void sync_deinit(void);

// Replace the communicator and refresh rank/size.
void sync_reinit(MPI_Comm comm);

// In-place global minimum.
void sync_min(void *buf, int num, MPI_Datatype type);

// In-place global maximum.
void sync_max(void *buf, int num, MPI_Datatype type);

// In-place global sum.
void sync_sum(void *buf, int num, MPI_Datatype type);

// Global dot product of two double arrays.
double sync_dot(const double *lhs, const double *rhs, int num);

// Global L2 norm of a double array.
double sync_norm(const double *buf, int num);

// In-place exclusive prefix sum; rank 0 receives zeros.
void sync_prefix(void *buf, int num, MPI_Datatype type);

// Build rank offsets with first entry zero and last entry equal to the global total.
void sync_offsets(const void *val, void *buf, int num, MPI_Datatype type);

// Gather values per rank into buffer.
void sync_gather(const void *val, void *buf, int num, MPI_Datatype type, int len);

// Rotate a fixed-capacity buffer and count to the next rank, receive from the previous rank.
void sync_rotate(void *buf, int *num, int cap, MPI_Datatype type, int len);

// Send to `dst` and receive from `src` in a single exchange.
void sync_exchange(const void *send, void *recv, int num_send, int num_recv, int dst, int src,
                   MPI_Datatype type, int len);

// Collect values at global indices from concatenated rank-local arrays.
void sync_collect(const void *send, void *recv, const long *global, int num_send, int num_recv,
                  MPI_Datatype type, int len);

// Post a non-blocking receive.
void sync_irecv(MPI_Request *request, void *buf, int num, int rank, int tag, MPI_Datatype type,
                int len);

// Post a non-blocking send.
void sync_isend(MPI_Request *request, const void *buf, int num, int rank, int tag,
                MPI_Datatype type, int len);

// Resize a derived datatype.
MPI_Datatype sync_resized(MPI_Datatype type, MPI_Aint extent);
