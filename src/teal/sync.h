#pragma once

#include <mpi.h>

#include "teal.h"

#define MPI_SCALAR (sizeof(scalar) == sizeof(float) ? MPI_FLOAT : MPI_DOUBLE)

enum { MPI_TAG_UB_MIN = 32767 };
#define sync_tag() ((__LINE__ % MPI_TAG_UB_MIN) + 1)

typedef struct {
    MPI_Comm comm;
    long rank;
    long size;
} Sync;

extern Sync sync;

void sync_init(int *argc, char ***argv);
void sync_reinit(MPI_Comm comm);
void sync_deinit(void);

long sync_lmin(long val);
long sync_lmax(long val);
long sync_lsum(long val);

scalar sync_fmin(scalar val);
scalar sync_fmax(scalar val);
scalar sync_fsum(scalar val);

vector sync_vector_min(vector val);
vector sync_vector_max(vector val);
vector sync_vector_sum(vector val);

long sync_exsum(long val);

scalar sync_dot(const scalar *lhs, const scalar *rhs, long num);
scalar sync_norm(const scalar *arr, long num);

MPI_Request *sync_irecv(const long *rank, const long *off, void *arr_, long num, long stride,
                        MPI_Datatype type, long tag);
MPI_Request *sync_isend(const long *rank, const long *off, const long *idx, const void *arr_,
                        long num, long stride, MPI_Datatype type, long tag);
