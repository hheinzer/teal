#pragma once

#include <mpi.h>

extern struct Teal {
    MPI_Comm comm;
    int rank, size;
} teal;

void teal_initialize(int *argc, char ***argv);

void teal_finalize(void);
