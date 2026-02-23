#include <assert.h>
#include <limits.h>
#include <stddef.h>
#include <string.h>

#include "sync2.h"
#include "teal2.h"
#include "utils2.h"

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

void sync2_collect(const void *send_, void *recv_, const long *idx_recv_, int num_send_,
                   int num_recv_, MPI_Datatype type)
{
    assert((send_ || num_send_ == 0) && num_send_ >= 0);
    assert(((recv_ && idx_recv_) || num_recv_ == 0) && num_recv_ >= 0);

    long *offset = teal2_calloc(sync2.size + 1, sizeof(*offset));
    MPI_Allgather(&(long){num_send_}, 1, MPI_LONG, &offset[1], 1, MPI_LONG, sync2.comm);
    for (int i = 0; i < sync2.size; i++) {
        offset[i + 1] += offset[i];
    }

    int *rank = teal2_calloc(num_recv_, sizeof(*rank));
    for (int i = 0; i < num_recv_; i++) {
        rank[i] = digitize(&idx_recv_[i], offset, sync2.size, sizeof(*offset), compare_long);
        assert(0 <= rank[i] && rank[i] < sync2.size);
    }

    int *num_recv = teal2_calloc(sync2.size, sizeof(*num_recv));
    for (int i = 0; i < num_recv_; i++) {
        num_recv[rank[i]] += 1;
    }

    int *off_recv = teal2_calloc(sync2.size + 1, sizeof(*off_recv));
    for (int i = 0; i < sync2.size; i++) {
        off_recv[i + 1] = off_recv[i] + num_recv[i];
    }

    int *cur_recv = teal2_calloc(sync2.size, sizeof(*cur_recv));
    copy(cur_recv, off_recv, sync2.size, sizeof(*cur_recv));

    int *idx_recv = teal2_calloc(num_recv_, sizeof(*idx_recv));
    for (int i = 0; i < num_recv_; i++) {
        long idx_local = idx_recv_[i] - offset[rank[i]];
        assert(idx_local <= INT_MAX);
        idx_recv[cur_recv[rank[i]]++] = (int)idx_local;
    }
    for (int i = 0; i < sync2.size; i++) {
        assert(cur_recv[i] == off_recv[i + 1]);
    }

    int *num_send = teal2_calloc(sync2.size, sizeof(*num_send));
    MPI_Alltoall(num_recv, 1, MPI_INT, num_send, 1, MPI_INT, sync2.comm);

    int *off_send = teal2_calloc(sync2.size + 1, sizeof(*off_send));
    for (int i = 0; i < sync2.size; i++) {
        off_send[i + 1] = off_send[i] + num_send[i];
    }

    int tot_send = off_send[sync2.size];
    int *idx_send = teal2_calloc(tot_send, sizeof(*idx_send));
    MPI_Alltoallv(idx_recv, num_recv, off_recv, MPI_INT, idx_send, num_send, off_send, MPI_INT,
                  sync2.comm);

    MPI_Aint lower_bound;
    MPI_Aint extent;
    MPI_Type_get_extent(type, &lower_bound, &extent);
    assert(lower_bound == 0 && 0 < extent && extent <= INT_MAX);

    char (*send)[extent] = teal2_calloc(tot_send, (int)extent);
    for (int i = 0; i < tot_send; i++) {
        assert(0 <= idx_send[i] && idx_send[i] < num_send_);
        memcpy(send[i], ((char (*)[extent])send_)[idx_send[i]], extent);
    }

    char (*recv)[extent] = teal2_calloc(num_recv_, (int)extent);
    MPI_Alltoallv(send, num_send, off_send, type, recv, num_recv, off_recv, type, sync2.comm);

    copy(cur_recv, off_recv, sync2.size, sizeof(*cur_recv));
    for (int i = 0; i < num_recv_; i++) {
        memcpy(((char (*)[extent])recv_)[i], recv[cur_recv[rank[i]]++], extent);
    }
    for (int i = 0; i < sync2.size; i++) {
        assert(cur_recv[i] == off_recv[i + 1]);
    }

    teal2_free(offset);
    teal2_free(rank);
    teal2_free(num_recv);
    teal2_free(off_recv);
    teal2_free(cur_recv);
    teal2_free(idx_recv);
    teal2_free(num_send);
    teal2_free(off_send);
    teal2_free(idx_send);
    teal2_free(send);
    teal2_free(recv);
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
