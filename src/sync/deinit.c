#include "sync.h"

void sync_deinit(void)
{
    MPI_Comm_free(&sync.comm);
    MPI_Finalize();
}
