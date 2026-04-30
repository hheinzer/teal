#include "sync.h"

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <string.h>

#include "teal.h"
#include "utils.h"

struct sync sync = {0};

void sync_init(int *argc, char ***argv)
{
    assert(argc && argv);
    MPI_Init(argc, argv);
    MPI_Comm_dup(MPI_COMM_WORLD, &sync.comm);
    MPI_Comm_rank(sync.comm, &sync.rank);
    MPI_Comm_size(sync.comm, &sync.size);
}

void sync_deinit(void)
{
    int flag;
    MPI_Initialized(&flag);
    if (flag) {
        MPI_Comm_free(&sync.comm);
        MPI_Finalize();
    }
}

void sync_reinit(MPI_Comm comm)
{
    int rank = sync.rank;

    MPI_Comm_free(&sync.comm);
    sync.comm = comm;

    MPI_Comm_rank(sync.comm, &sync.rank);
    MPI_Comm_size(sync.comm, &sync.size);

    if (rank != sync.rank) {
        teal_verbose("rank remapped (%d -> %d)", rank, sync.rank);
    }
}

void sync_min(void *buf, int num, MPI_Datatype type)
{
    assert(buf && num > 0);
    MPI_Allreduce(MPI_IN_PLACE, buf, num, type, MPI_MIN, sync.comm);
}

void sync_max(void *buf, int num, MPI_Datatype type)
{
    assert(buf && num > 0);
    MPI_Allreduce(MPI_IN_PLACE, buf, num, type, MPI_MAX, sync.comm);
}

void sync_sum(void *buf, int num, MPI_Datatype type)
{
    assert(buf && num > 0);
    MPI_Allreduce(MPI_IN_PLACE, buf, num, type, MPI_SUM, sync.comm);
}

double sync_dot(const double *lhs, const double *rhs, int num)
{
    double dot = 0;
    for (int i = 0; i < num; i++) {
        dot += lhs[i] * rhs[i];
    }
    sync_sum(&dot, 1, MPI_DOUBLE);
    return dot;
}

double sync_norm(const double *buf, int num)
{
    return sqrt(sync_dot(buf, buf, num));
}

void sync_prefix(void *buf, int num, MPI_Datatype type)
{
    assert(buf && num > 0);

    MPI_Exscan(MPI_IN_PLACE, buf, num, type, MPI_SUM, sync.comm);

    if (sync.rank == 0) {
        int size;
        MPI_Type_size(type, &size);
        assert(size > 0);

        memset(buf, 0, (size_t)num * size);
    }
}

void sync_offsets(const void *val, void *buf_, int num, MPI_Datatype type)
{
    assert(val && buf_ && num > 0);

    int size;
    MPI_Type_size(type, &size);
    assert(size > 0);

    char (*buf)[num * size] = buf_;
    memset(buf[0], 0, sizeof(*buf));

    MPI_Scan(val, buf[sync.rank + 1], num, type, MPI_SUM, sync.comm);
    MPI_Allgather(MPI_IN_PLACE, num, type, buf[1], num, type, sync.comm);
}

void sync_gather(const void *val, void *buf, int num, MPI_Datatype type, int len)
{
    assert(buf && num > 0 && len > 0);
    MPI_Allgather(val, num * len, type, buf, num * len, type, sync.comm);
}

void sync_rotate(void *buf, int *num, int cap, MPI_Datatype type, int len)
{
    assert(buf && num && cap > 0 && len > 0);
    int next = (sync.rank + 1) % sync.size;
    int prev = (sync.rank - 1 + sync.size) % sync.size;
    MPI_Sendrecv_replace(num, 1, MPI_INT, next, 0, prev, 0, sync.comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv_replace(buf, cap * len, type, next, 1, prev, 1, sync.comm, MPI_STATUS_IGNORE);
}

void sync_exchange(const void *send, void *recv, int num_send, int num_recv, int dst, int src,
                   MPI_Datatype type, int len)
{
    assert(send && recv && num_send > 0 && num_recv > 0 && len > 0);
    MPI_Sendrecv(send, num_send * len, type, dst, 0, recv, num_recv * len, type, src, 0, sync.comm,
                 MPI_STATUS_IGNORE);
}

void sync_collect(const void *send_, void *recv_, const long *global, int num_send_, int num_recv_,
                  MPI_Datatype type, int len)
{
    assert((send_ || num_send_ == 0) && num_send_ >= 0);
    assert(((recv_ && global) || num_recv_ == 0) && num_recv_ >= 0);
    assert(len > 0);

    long *offset = teal_calloc(sync.size + 1, sizeof(*offset));
    MPI_Allgather(&(long){num_send_}, 1, MPI_LONG, &offset[1], 1, MPI_LONG, sync.comm);
    for (int i = 0; i < sync.size; i++) {
        offset[i + 1] += offset[i];
    }

    int *rank = teal_calloc(num_recv_, sizeof(*rank));
    for (int i = 0; i < num_recv_; i++) {
        rank[i] = digitize(&global[i], offset, sync.size, sizeof(*offset), compare_long);
        assert(0 <= rank[i] && rank[i] < sync.size);
    }

    int *num_recv = teal_calloc(sync.size, sizeof(*num_recv));
    for (int i = 0; i < num_recv_; i++) {
        num_recv[rank[i]] += 1;
    }

    int *off_recv = teal_calloc(sync.size + 1, sizeof(*off_recv));
    for (int i = 0; i < sync.size; i++) {
        off_recv[i + 1] = off_recv[i] + num_recv[i];
    }

    int *cur_recv = teal_calloc(sync.size, sizeof(*cur_recv));
    copy(cur_recv, off_recv, sync.size, sizeof(*cur_recv));

    int *idx_recv = teal_calloc(num_recv_, sizeof(*idx_recv));
    for (int i = 0; i < num_recv_; i++) {
        long local = global[i] - offset[rank[i]];
        assert(local <= INT_MAX);
        idx_recv[cur_recv[rank[i]]++] = (int)local;
    }
    for (int i = 0; i < sync.size; i++) {
        assert(cur_recv[i] == off_recv[i + 1]);
    }

    int *num_send = teal_calloc(sync.size, sizeof(*num_send));
    MPI_Alltoall(num_recv, 1, MPI_INT, num_send, 1, MPI_INT, sync.comm);

    int *off_send = teal_calloc(sync.size + 1, sizeof(*off_send));
    for (int i = 0; i < sync.size; i++) {
        off_send[i + 1] = off_send[i] + num_send[i];
    }

    int tot_send = off_send[sync.size];
    int *idx_send = teal_calloc(tot_send, sizeof(*idx_send));
    MPI_Alltoallv(idx_recv, num_recv, off_recv, MPI_INT, idx_send, num_send, off_send, MPI_INT,
                  sync.comm);

    if (len > 1) {
        MPI_Type_contiguous(len, type, &type);
        MPI_Type_commit(&type);
    }

    int size;
    MPI_Type_size(type, &size);
    assert(size > 0);

    char (*send)[size] = teal_calloc(tot_send, size);
    for (int i = 0; i < tot_send; i++) {
        assert(0 <= idx_send[i] && idx_send[i] < num_send_);
        memcpy(send[i], ((char (*)[size])send_)[idx_send[i]], size);
    }

    char (*recv)[size] = teal_calloc(num_recv_, size);
    MPI_Alltoallv(send, num_send, off_send, type, recv, num_recv, off_recv, type, sync.comm);

    if (len > 1) {
        MPI_Type_free(&type);
    }

    copy(cur_recv, off_recv, sync.size, sizeof(*cur_recv));
    for (int i = 0; i < num_recv_; i++) {
        memcpy(((char (*)[size])recv_)[i], recv[cur_recv[rank[i]]++], size);
    }
    for (int i = 0; i < sync.size; i++) {
        assert(cur_recv[i] == off_recv[i + 1]);
    }

    teal_free(offset);
    teal_free(rank);
    teal_free(num_recv);
    teal_free(off_recv);
    teal_free(cur_recv);
    teal_free(idx_recv);
    teal_free(num_send);
    teal_free(off_send);
    teal_free(idx_send);
    teal_free(send);
    teal_free(recv);
}

void sync_irecv(MPI_Request *request, void *buf, int num, int rank, int tag, MPI_Datatype type,
                int len)
{
    assert(buf && num > 0 && len > 0);
    MPI_Irecv(buf, num * len, type, rank, tag, sync.comm, request);
}

void sync_isend(MPI_Request *request, const void *buf, int num, int rank, int tag,
                MPI_Datatype type, int len)
{
    assert(buf && num > 0 && len > 0);
    MPI_Isend(buf, num * len, type, rank, tag, sync.comm, request);
}

MPI_Datatype sync_resized(MPI_Datatype type, MPI_Aint extent)
{
    assert(extent > 0);
    MPI_Datatype resized;
    MPI_Type_create_resized(type, 0, extent, &resized);
    MPI_Type_commit(&resized);
    MPI_Type_free(&type);
    return resized;
}
