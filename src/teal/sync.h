#pragma once

#include <mpi.h>

#include "teal.h"

typedef struct {
    MPI_Comm comm;
    int rank;
    int size;
} Sync;

extern Sync sync;

extern MPI_Datatype vector_type;

void sync_init(int *argc, char ***argv);
void sync_reinit(MPI_Comm comm);

int sync_tag(void);

long sync_lmin(long val);
long sync_lmax(long val);
long sync_lsum(long val);

double sync_fmin(double val);
double sync_fmax(double val);
double sync_fsum(double val);

vector sync_vmin(vector val);
vector sync_vmax(vector val);
vector sync_vsum(vector val);

long sync_lexsum(long val);
double sync_fexsum(double val);

void sync_finalize(void);
