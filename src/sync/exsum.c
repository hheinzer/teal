#include "sync.h"

long sync_exsum(long val)
{
    long exsum = 0;
    MPI_Exscan(&val, &exsum, 1, MPI_LONG, MPI_SUM, sync.comm);
    return (sync.rank == 0) ? 0 : exsum;
}
