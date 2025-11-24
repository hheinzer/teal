#include <mpi.h>
#include <stdlib.h>

#include "teal.h"

void teal_exit(int status)
{
    int flag;
    MPI_Initialized(&flag);
    if (flag) {
        MPI_Finalize();
    }
    exit(status);
}
