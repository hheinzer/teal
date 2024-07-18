#include "sync.h"

#include <mpi.h>
#include <stdlib.h>

#include "memory.h"
#include "mesh.h"
#include "teal.h"
#include "utils.h"

static void irecv(const int size, const long *j_recv, const long nu, double (*u)[nu],
                  MPI_Request *req);
static void isend(const int size, const long *i_send, const long *send, const long nu,
                  const double (*u)[nu], double (*buf)[nu], MPI_Request *req);
static void wait(const int size, MPI_Request *req);

long sync_max_long(const long v)
{
    long r = 0;
    MPI_Allreduce(&v, &r, 1, MPI_LONG, MPI_MAX, teal.comm);
    return r;
}
double sync_max_double(const double v)
{
    double r = 0;
    MPI_Allreduce(&v, &r, 1, MPI_DOUBLE, MPI_MAX, teal.comm);
    return r;
}

long sync_min_long(const long v)
{
    long r = 0;
    MPI_Allreduce(&v, &r, 1, MPI_LONG, MPI_MIN, teal.comm);
    return r;
}
double sync_min_double(const double v)
{
    double r = 0;
    MPI_Allreduce(&v, &r, 1, MPI_DOUBLE, MPI_MIN, teal.comm);
    return r;
}

long sync_sum_long(const long v)
{
    long r = 0;
    MPI_Allreduce(&v, &r, 1, MPI_LONG, MPI_SUM, teal.comm);
    return r;
}
double sync_sum_double(const double v)
{
    double r = 0;
    MPI_Allreduce(&v, &r, 1, MPI_DOUBLE, MPI_SUM, teal.comm);
    return r;
}

long sync_exscan_sum_long(const long v)
{
    long r = 0;
    MPI_Exscan(&v, &r, 1, MPI_LONG, MPI_SUM, teal.comm);
    return r;
}
double sync_exscan_sum_double(const double v)
{
    double r = 0;
    MPI_Exscan(&v, &r, 1, MPI_DOUBLE, MPI_SUM, teal.comm);
    return r;
}

long *sync_gather_long(long *v, const long n)
{
    cleanup int *counts = memory_calloc(teal.size, sizeof(*counts));
    cleanup int *displs = memory_calloc(teal.size, sizeof(*displs));
    counts[teal.rank] = n;
    displs[teal.rank] = sync_exscan_sum(n);
    ensure(counts[teal.rank] > 0 && displs[teal.rank] > 0);

    MPI_Allgather(&counts[teal.rank], 1, MPI_INT, counts, 1, MPI_INT, teal.comm);
    MPI_Allgather(&displs[teal.rank], 1, MPI_INT, displs, 1, MPI_INT, teal.comm);

    long *buf = memory_calloc(sync_sum(n), sizeof(*buf));
    MPI_Allgatherv(v, n, MPI_LONG, buf, counts, displs, MPI_LONG, teal.comm);
    free(v);

    return buf;
}
double *sync_gather_double(double *v, const long n)
{
    cleanup int *counts = memory_calloc(teal.size, sizeof(*counts));
    cleanup int *displs = memory_calloc(teal.size, sizeof(*displs));
    counts[teal.rank] = n;
    displs[teal.rank] = sync_exscan_sum(n);
    ensure(counts[teal.rank] >= 0 && displs[teal.rank] >= 0);

    MPI_Allgather(&counts[teal.rank], 1, MPI_INT, counts, 1, MPI_INT, teal.comm);
    MPI_Allgather(&displs[teal.rank], 1, MPI_INT, displs, 1, MPI_INT, teal.comm);

    double *buf = memory_calloc(sync_sum(n), sizeof(*buf));
    MPI_Allgatherv(v, n, MPI_DOUBLE, buf, counts, displs, MPI_DOUBLE, teal.comm);
    free(v);

    return buf;
}

void sync_all(const Mesh *mesh, const long nu, double (*u)[nu])
{
    const ALIAS(j_recv, mesh->sync.j_recv);
    const ALIAS(i_send, mesh->sync.i_send);
    const ALIAS(send, mesh->sync.send);

    cleanup double(*buf)[nu] = memory_calloc(i_send[teal.size], sizeof(*buf));
    cleanup MPI_Request(*req)[teal.size] = memory_calloc(2, sizeof(*req));

    irecv(teal.size, j_recv, nu, u, req[0]);
    isend(teal.size, i_send, send, nu, u, buf, req[1]);
    wait(2 * teal.size, *req);
}

void sync_begin(Equations *eqns, const long nu, double (*u)[nu])
{
    const ALIAS(j_recv, eqns->mesh->sync.j_recv);
    irecv(teal.size, j_recv, nu, u, eqns->sync.recv);

    const ALIAS(i_send, eqns->mesh->sync.i_send);
    const ALIAS(send, eqns->mesh->sync.send);
    double(*buf)[nu] = (void *)eqns->sync.buf;
    isend(teal.size, i_send, send, nu, u, buf, eqns->sync.send);
}

void sync_wait(Equations *eqns)
{
    wait(teal.size, eqns->sync.recv);
}

void sync_end(Equations *eqns)
{
    wait(teal.size, eqns->sync.send);
}

static void irecv(const int size, const long *j_recv, const long nu, double (*u)[nu],
                  MPI_Request *req)
{
    for (long rank = 0; rank < size; ++rank) {
        const long count = (j_recv[rank + 1] - j_recv[rank]) * nu;
        MPI_Irecv(&u[j_recv[rank]], count, MPI_DOUBLE, rank, 0, teal.comm, &req[rank]);
    }
}

static void isend(const int size, const long *i_send, const long *send, const long nu,
                  const double (*u)[nu], double (*buf)[nu], MPI_Request *req)
{
    for (long rank = 0; rank < size; ++rank)
        for (long i = i_send[rank]; i < i_send[rank + 1]; ++i)
            for (long j = 0; j < nu; ++j) buf[i][j] = u[send[i]][j];

    for (long rank = 0; rank < size; ++rank) {
        const long count = (i_send[rank + 1] - i_send[rank]) * nu;
        MPI_Isend(buf[i_send[rank]], count, MPI_DOUBLE, rank, 0, teal.comm, &req[rank]);
    }
}

static void wait(const int size, MPI_Request *req)
{
    MPI_Waitall(size, req, MPI_STATUSES_IGNORE);
}
