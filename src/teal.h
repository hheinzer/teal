#pragma once

#include <mpi.h>

extern struct Teal {
    MPI_Comm comm;
    int rank, size;
    int quiet;
    const char *restart;
} teal;

/* Initialize teal and print startup message. This must be called before any call to the other
 * functions of teal. */
void teal_initialize(int *argc, char ***argv);

/* Finalize teal and print end message. This must be called when you are done using teal. */
void teal_finalize(void);
