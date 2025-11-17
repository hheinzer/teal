#include "sync.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "assert.h"
#include "teal/arena.h"
#include "teal/array.h"

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

number sync_lmin(number val)
{
    number min = val;
    MPI_Allreduce(&val, &min, 1, MPI_NUMBER, MPI_MIN, sync.comm);
    return min;
}

number sync_lmax(number val)
{
    number max = val;
    MPI_Allreduce(&val, &max, 1, MPI_NUMBER, MPI_MAX, sync.comm);
    return max;
}

number sync_lsum(number val)
{
    number sum = 0;
    MPI_Allreduce(&val, &sum, 1, MPI_NUMBER, MPI_SUM, sync.comm);
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

number sync_lexsum(number val)
{
    number exsum = 0;
    MPI_Exscan(&val, &exsum, 1, MPI_NUMBER, MPI_SUM, sync.comm);
    return (sync.rank == 0) ? 0 : exsum;
}

scalar sync_fdot(const scalar *lhs, const scalar *rhs, number num)
{
    return sync_fsum(array_fdot(lhs, rhs, num));
}

scalar sync_fnorm(const scalar *arr, number num)
{
    return sqrt(sync_fdot(arr, arr, num));
}

MPI_Request *sync_irecv_scalar(const number *rank, const number *off, void *arr_, number num,
                               number stride, int tag)
{
    scalar(*arr)[stride] = arr_;
    MPI_Request *req = arena_malloc(num, sizeof(*req));
    for (number i = 0; i < num; i++) {
        number count = (off[i + 1] - off[i]) * stride;
        MPI_Irecv(arr[off[i]], count, MPI_SCALAR, rank[i], tag, sync.comm, &req[i]);
    }
    return req;
}

MPI_Request *sync_isend_scalar(const number *rank, const number *off, const number *idx,
                               const void *arr_, number num, number stride, int tag)
{
    const scalar(*arr)[stride] = arr_;
    scalar(*buf)[stride] = arena_malloc(off[num], sizeof(*buf));
    for (number i = 0; i < num; i++) {
        for (number j = off[i]; j < off[i + 1]; j++) {
            memcpy(buf[j], arr[idx[j]], sizeof(*arr));
        }
    }
    MPI_Request *req = arena_malloc(num, sizeof(*req));
    for (number i = 0; i < num; i++) {
        number count = (off[i + 1] - off[i]) * stride;
        MPI_Isend(buf[off[i]], count, MPI_SCALAR, rank[i], tag, sync.comm, &req[i]);
    }
    return req;
}

MPI_Request *sync_irecv_vector(const number *rank, const number *off, void *arr_, number num,
                               number stride, int tag)
{
    return sync_irecv_scalar(rank, off, arr_, num, stride * 3, tag);
}

MPI_Request *sync_isend_vector(const number *rank, const number *off, const number *idx,
                               const void *arr_, number num, number stride, int tag)
{
    return sync_isend_scalar(rank, off, idx, arr_, num, stride * 3, tag);
}

void sync_finalize(void)
{
    MPI_Comm_free(&sync.comm);
    MPI_Finalize();
}
