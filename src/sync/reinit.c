#include <assert.h>

#include "sync.h"

void sync_reinit(MPI_Comm comm)
{
    assert(comm != MPI_COMM_NULL);
    MPI_Comm_free(&sync.comm);
    sync.comm = comm;
    MPI_Comm_rank(sync.comm, &sync.rank);
    MPI_Comm_size(sync.comm, &sync.size);
}
