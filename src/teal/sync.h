#pragma once

#include <mpi.h>
#include <stdint.h>

#include "teal.h"

#define MPI_NUMBER (sizeof(number) == sizeof(int32_t) ? MPI_INT32_T : MPI_INT64_T)
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

number sync_lmin(number val);
number sync_lmax(number val);
number sync_lsum(number val);

scalar sync_fmin(scalar val);
scalar sync_fmax(scalar val);
scalar sync_fsum(scalar val);

vector sync_vmin(vector val);
vector sync_vmax(vector val);
vector sync_vsum(vector val);

/* Return exclusive prefix sum of `val`. */
number sync_lexsum(number val);

void sync_finalize(void);
