#include "teal.h"

#include <stdio.h>
#include <time.h>

#include "core/utils.h"

struct Teal teal = {0};

void teal_initialize(int *argc, char ***argv)
{
    MPI_Init(argc, argv);

    teal.comm = MPI_COMM_WORLD;
    MPI_Comm_rank(teal.comm, &teal.rank);
    MPI_Comm_size(teal.comm, &teal.size);

    char buf[128];
    strftime(buf, sizeof(buf), "%a %b %e %T %Y", localtime((time_t[]){time(0)}));

    if (teal.rank == 0) {
        printf("Hello, World! This is teal!\n");
        printf(" | " KEYFMT ": %s\n", "start time", buf);
        printf(" | " KEYFMT ": %d\n", "number of ranks", teal.size);
    }
}

void teal_finalize(void)
{
    char buf[128];
    strftime(buf, sizeof(buf), "%a %b %e %T %Y", localtime((time_t[]){time(0)}));

    if (teal.rank == 0) {
        printf("Goodbye, World!\n");
        printf(" | " KEYFMT ": %s\n", "end time", buf);
    }

    MPI_Finalize();
}
