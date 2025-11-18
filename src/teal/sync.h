#pragma once

#include <mpi.h>
#include <stdint.h>

#include "teal.h"

#define MPI_SCALAR (sizeof(scalar) == sizeof(float) ? MPI_FLOAT : MPI_DOUBLE)

enum { MPI_TAG_UB_MIN = 32767 };
#define sync_tag() ((__LINE__ % MPI_TAG_UB_MIN) + 1)

typedef struct {
    MPI_Comm comm;
    int rank;
    int size;
} Sync;

extern Sync sync;

void sync_init(int *argc, char ***argv);
void sync_reinit(MPI_Comm comm);

void sync_exit(int status) __attribute((noreturn));
void sync_abort(void) __attribute((noreturn));

int sync_lmin(int val);
int sync_lmax(int val);
int sync_lsum(int val);

scalar sync_fmin(scalar val);
scalar sync_fmax(scalar val);
scalar sync_fsum(scalar val);

vector sync_vector_min(vector val);
vector sync_vector_max(vector val);
vector sync_vector_sum(vector val);

int sync_lexsum(int val);

scalar sync_fdot(const scalar *lhs, const scalar *rhs, int num);
scalar sync_fnorm(const scalar *arr, int num);

MPI_Request *sync_irecv(const int *rank, const int *off, void *arr_, int num, int stride,
                        MPI_Datatype type, int tag);
MPI_Request *sync_isend(const int *rank, const int *off, const int *idx, const void *arr_, int num,
                        int stride, MPI_Datatype type, int tag);

void sync_finalize(void);
