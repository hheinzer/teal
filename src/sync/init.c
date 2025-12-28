#include <assert.h>

#include "sync.h"

struct Sync sync = {0};

void sync_init(int *argc, char ***argv)
{
    assert(argc && argv);
    MPI_Init(argc, argv);
    MPI_Comm_dup(MPI_COMM_WORLD, &sync.comm);
    MPI_Comm_rank(sync.comm, &sync.rank);
    MPI_Comm_size(sync.comm, &sync.size);
}
