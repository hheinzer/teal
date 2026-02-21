#include <assert.h>
#include <stddef.h>
#include <string.h>

#include "sync2.h"
#include "teal2.h"

struct sync sync2 = {0};

void sync2_init(int *argc, char ***argv)
{
    assert(argc && argv);
    MPI_Init(argc, argv);
    MPI_Comm_dup(MPI_COMM_WORLD, &sync2.comm);
    MPI_Comm_rank(sync2.comm, &sync2.rank);
    MPI_Comm_size(sync2.comm, &sync2.size);
}

void sync2_deinit(void)
{
    int flag;
    MPI_Initialized(&flag);
    if (flag) {
        MPI_Comm_free(&sync2.comm);
        MPI_Finalize();
    }
}

void sync2_reinit(MPI_Comm comm)
{
    int rank = sync2.rank;

    MPI_Comm_free(&sync2.comm);
    sync2.comm = comm;

    MPI_Comm_rank(sync2.comm, &sync2.rank);
    MPI_Comm_size(sync2.comm, &sync2.size);

    if (rank != sync2.rank) {
        teal2_verbose("rank remapped (%d -> %d)", rank, sync2.rank);
    }
}

void sync2_min(void *buf, int num, MPI_Datatype type)
{
    assert(buf && num > 0);
    MPI_Allreduce(MPI_IN_PLACE, buf, num, type, MPI_MIN, sync2.comm);
}

void sync2_max(void *buf, int num, MPI_Datatype type)
{
    assert(buf && num > 0);
    MPI_Allreduce(MPI_IN_PLACE, buf, num, type, MPI_MAX, sync2.comm);
}

void sync2_sum(void *buf, int num, MPI_Datatype type)
{
    assert(buf && num > 0);
    MPI_Allreduce(MPI_IN_PLACE, buf, num, type, MPI_SUM, sync2.comm);
}

void sync2_prefix(void *buf, int num, MPI_Datatype type)
{
    assert(buf && num > 0);
    MPI_Exscan(MPI_IN_PLACE, buf, num, type, MPI_SUM, sync2.comm);
    if (sync2.rank == 0) {
        int size;
        MPI_Type_size(type, &size);
        assert(size > 0);
        memset(buf, 0, (size_t)num * size);
    }
}

void sync2_gather(const void *val, void *buf, int num, MPI_Datatype type)
{
    assert(buf && num > 0);
    MPI_Allgather(val, num, type, buf, num, type, sync2.comm);
}

void sync2_offsets(const void *val, void *buf_, int num, MPI_Datatype type)
{
    assert(val && buf_ && num > 0);
    int size;
    MPI_Type_size(type, &size);
    assert(size > 0);
    char (*buf)[num * size] = buf_;
    memset(buf[0], 0, sizeof(*buf));
    MPI_Scan(val, buf[sync2.rank + 1], num, type, MPI_SUM, sync2.comm);
    MPI_Allgather(MPI_IN_PLACE, num, type, buf[1], num, type, sync2.comm);
}

void sync2_rotate(void *buf, int *num, int cap, MPI_Datatype type)
{
    assert(buf && num && cap > 0);
    int next = (sync2.rank + 1) % sync2.size;
    int prev = (sync2.rank - 1 + sync2.size) % sync2.size;
    MPI_Sendrecv_replace(num, 1, MPI_INT, next, 0, prev, 0, sync2.comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv_replace(buf, cap, type, next, 1, prev, 1, sync2.comm, MPI_STATUS_IGNORE);
}

void sync2_exchange(const void *send, void *recv, int num_send, int num_recv, int dst, int src,
                    MPI_Datatype type)
{
    assert(send && recv && num_send > 0 && num_recv > 0);
    MPI_Sendrecv(send, num_send, type, dst, 0, recv, num_recv, type, src, 0, sync2.comm,
                 MPI_STATUS_IGNORE);
}

MPI_Datatype sync2_resized(MPI_Datatype type, MPI_Aint extent)
{
    assert(extent > 0);
    MPI_Datatype resized;
    MPI_Type_create_resized(type, 0, extent, &resized);
    MPI_Type_commit(&resized);
    MPI_Type_free(&type);
    return resized;
}
