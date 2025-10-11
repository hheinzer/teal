#pragma once

#include <mpi.h>

#include "teal.h"

#define MPI_SCALAR (sizeof(scalar) == sizeof(float) ? MPI_FLOAT : MPI_DOUBLE)

typedef struct {
    MPI_Comm comm;
    int rank;
    int size;
} Sync;

extern Sync sync;

void sync_init(int *argc, char ***argv);
void sync_reinit(MPI_Comm comm);

int sync_tag(void);

long sync_lmin(long val);
long sync_lmax(long val);
long sync_lsum(long val);

scalar sync_fmin(scalar val);
scalar sync_fmax(scalar val);
scalar sync_fsum(scalar val);

vector sync_vmin(vector val);
vector sync_vmax(vector val);
vector sync_vsum(vector val);

/* Return exclusive prefix sum of `val`. */
long sync_lexsum(long val);

void sync_finalize(void);
