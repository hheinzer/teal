#include "sync.h"

#include <stdlib.h>

#include "assert.h"

Sync sync = {0};

void sync_init(int *argc, char ***argv)
{
    assert(argc && argv);
    MPI_Init(argc, argv);
    MPI_Comm_dup(MPI_COMM_WORLD, &sync.comm);
    MPI_Comm_rank(sync.comm, &sync.rank);
    MPI_Comm_size(sync.comm, &sync.size);
}

void sync_reinit(MPI_Comm comm)
{
    assert(comm != MPI_COMM_NULL);
    MPI_Comm_free(&sync.comm);
    sync.comm = comm;
    MPI_Comm_rank(sync.comm, &sync.rank);
    MPI_Comm_size(sync.comm, &sync.size);
}

void sync_exit(int status)
{
    int flag;
    MPI_Initialized(&flag);
    if (flag) {
        MPI_Finalize();
    }
    exit(status);
}

void sync_abort(void)
{
    int flag;
    MPI_Initialized(&flag);
    if (flag) {
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    abort();
}

long sync_lmin(long val)
{
    long min = val;
    MPI_Allreduce(&val, &min, 1, MPI_LONG, MPI_MIN, sync.comm);
    return min;
}

long sync_lmax(long val)
{
    long max = val;
    MPI_Allreduce(&val, &max, 1, MPI_LONG, MPI_MAX, sync.comm);
    return max;
}

long sync_lsum(long val)
{
    long sum = 0;
    MPI_Allreduce(&val, &sum, 1, MPI_LONG, MPI_SUM, sync.comm);
    return sum;
}

scalar sync_fmin(scalar val)
{
    scalar min = val;
    MPI_Allreduce(&val, &min, 1, MPI_SCALAR, MPI_MIN, sync.comm);
    return min;
}

scalar sync_fmax(scalar val)
{
    scalar max = val;
    MPI_Allreduce(&val, &max, 1, MPI_SCALAR, MPI_MAX, sync.comm);
    return max;
}

scalar sync_fsum(scalar val)
{
    scalar sum = 0;
    MPI_Allreduce(&val, &sum, 1, MPI_SCALAR, MPI_SUM, sync.comm);
    return sum;
}

vector sync_vmin(vector val)
{
    vector min = val;
    MPI_Allreduce(&val, &min, 3, MPI_SCALAR, MPI_MIN, sync.comm);
    return min;
}

vector sync_vmax(vector val)
{
    vector max = val;
    MPI_Allreduce(&val, &max, 3, MPI_SCALAR, MPI_MAX, sync.comm);
    return max;
}

vector sync_vsum(vector val)
{
    vector sum = {0};
    MPI_Allreduce(&val, &sum, 3, MPI_SCALAR, MPI_SUM, sync.comm);
    return sum;
}

long sync_lexsum(long val)
{
    long exsum = 0;
    MPI_Exscan(&val, &exsum, 1, MPI_LONG, MPI_SUM, sync.comm);
    return (sync.rank == 0) ? 0 : exsum;
}

void sync_finalize(void)
{
    MPI_Comm_free(&sync.comm);
    MPI_Finalize();
}
