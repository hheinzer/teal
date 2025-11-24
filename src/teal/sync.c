#include "sync.h"

#include <assert.h>
#include <math.h>
#include <string.h>

#include "teal/arena.h"
#include "teal/array.h"

Sync sync = {0};

void sync_init(int *argc, char ***argv)
{
    assert(argc && argv);

    MPI_Init(argc, argv);
    MPI_Comm_dup(MPI_COMM_WORLD, &sync.comm);

    int rank;
    MPI_Comm_rank(sync.comm, &rank);
    sync.rank = rank;

    int size;
    MPI_Comm_size(sync.comm, &size);
    sync.size = size;
}

void sync_reinit(MPI_Comm comm)
{
    assert(comm != MPI_COMM_NULL);

    MPI_Comm_free(&sync.comm);
    sync.comm = comm;

    int rank;
    MPI_Comm_rank(sync.comm, &rank);
    sync.rank = rank;

    int size;
    MPI_Comm_size(sync.comm, &size);
    sync.size = size;
}

void sync_deinit(void)
{
    MPI_Comm_free(&sync.comm);
    MPI_Finalize();
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

vector sync_vector_min(vector val)
{
    vector min = val;
    MPI_Allreduce(&val, &min, 3, MPI_SCALAR, MPI_MIN, sync.comm);
    return min;
}

vector sync_vector_max(vector val)
{
    vector max = val;
    MPI_Allreduce(&val, &max, 3, MPI_SCALAR, MPI_MAX, sync.comm);
    return max;
}

vector sync_vector_sum(vector val)
{
    vector sum = {0};
    MPI_Allreduce(&val, &sum, 3, MPI_SCALAR, MPI_SUM, sync.comm);
    return sum;
}

long sync_exsum(long val)
{
    long exsum = 0;
    MPI_Exscan(&val, &exsum, 1, MPI_LONG, MPI_SUM, sync.comm);
    return (sync.rank == 0) ? 0 : exsum;
}

scalar sync_dot(const scalar *lhs, const scalar *rhs, long num)
{
    return sync_fsum(array_dot(lhs, rhs, num));
}

scalar sync_norm(const scalar *arr, long num)
{
    return sqrt(sync_dot(arr, arr, num));
}

MPI_Request *sync_irecv(const long *rank, const long *off, void *arr_, long num, long stride,
                        MPI_Datatype type, long tag)
{
    int size;
    MPI_Type_size(type, &size);
    char (*arr)[stride * size] = arr_;
    MPI_Request *req = arena_malloc(num, sizeof(*req));
    for (long i = 0; i < num; i++) {
        long count = (off[i + 1] - off[i]) * stride;
        MPI_Irecv(arr[off[i]], count, type, rank[i], tag, sync.comm, &req[i]);
    }
    return req;
}

MPI_Request *sync_isend(const long *rank, const long *off, const long *idx, const void *arr_,
                        long num, long stride, MPI_Datatype type, long tag)
{
    int size;
    MPI_Type_size(type, &size);
    const char (*arr)[stride * size] = arr_;
    char (*buf)[stride * size] = arena_malloc(off[num], sizeof(*buf));
    for (long i = 0; i < num; i++) {
        for (long j = off[i]; j < off[i + 1]; j++) {
            memcpy(buf[j], arr[idx[j]], sizeof(*arr));
        }
    }
    MPI_Request *req = arena_malloc(num, sizeof(*req));
    for (long i = 0; i < num; i++) {
        long count = (off[i + 1] - off[i]) * stride;
        MPI_Isend(buf[off[i]], count, type, rank[i], tag, sync.comm, &req[i]);
    }
    return req;
}
