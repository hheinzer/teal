#include "teal.h"

#include <mpi.h>
#include <stdio.h>
#include <time.h>

#include "global.h"

void teal_init(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        char buf[128] = "";
        strftime(buf, sizeof(buf), "%a %b %e %T %Y", localtime((time_t[]){time(0)}));

        printf("Hello, World! This is teal!\n");
        printf(" | " FMT_KEY ": %s\n", "start time", buf);
        printf(" | " FMT_KEY ": %d\n", "mpi size", size);
    }
}

void teal_finalize(void)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        char buf[128] = "";
        strftime(buf, sizeof(buf), "%a %b %e %T %Y", localtime((time_t[]){time(0)}));

        printf("Goodbye, World!\n");
        printf(" | " FMT_KEY ": %s\n", "end time", buf);
    }

    MPI_Finalize();
}
