#include "sync.h"

#include <mpi.h>

#include "memory.h"
#include "mesh.h"
#include "teal.h"
#include "utils.h"

static void irecv(const long *j_recv, double *u, MPI_Request *req, long nu);

static void isend(const long *i_send, const long *send, const double *u, double *buf,
                  MPI_Request *req, long nu);

static void wait(MPI_Request *req);

long x__sync_max_long(long v)
{
    long r = 0;
    MPI_Allreduce(&v, &r, 1, MPI_LONG, MPI_MAX, teal.comm);
    return r;
}
double x__sync_max_double(double v)
{
    double r = 0;
    MPI_Allreduce(&v, &r, 1, MPI_DOUBLE, MPI_MAX, teal.comm);
    return r;
}

long x__sync_min_long(long v)
{
    long r = 0;
    MPI_Allreduce(&v, &r, 1, MPI_LONG, MPI_MIN, teal.comm);
    return r;
}
double x__sync_min_double(double v)
{
    double r = 0;
    MPI_Allreduce(&v, &r, 1, MPI_DOUBLE, MPI_MIN, teal.comm);
    return r;
}

long x__sync_sum_long(long v)
{
    long r = 0;
    MPI_Allreduce(&v, &r, 1, MPI_LONG, MPI_SUM, teal.comm);
    return r;
}
double x__sync_sum_double(double v)
{
    double r = 0;
    MPI_Allreduce(&v, &r, 1, MPI_DOUBLE, MPI_SUM, teal.comm);
    return r;
}

long x__sync_exscan_sum_long(long v)
{
    long r = 0;
    MPI_Exscan(&v, &r, 1, MPI_LONG, MPI_SUM, teal.comm);
    return r;
}
double x__sync_exscan_sum_double(double v)
{
    double r = 0;
    MPI_Exscan(&v, &r, 1, MPI_DOUBLE, MPI_SUM, teal.comm);
    return r;
}

void sync_all(const Mesh *mesh, double *u, long nu)
{
    const ALIAS(j_recv, mesh->sync.j_recv);
    const ALIAS(i_send, mesh->sync.i_send);
    const ALIAS(send, mesh->sync.send);

    smart double *buf = memory_calloc(i_send[teal.size] * nu, sizeof(*buf));
    smart MPI_Request *req_recv = memory_calloc(teal.size, sizeof(*req_recv));  // NOLINT
    smart MPI_Request *req_send = memory_calloc(teal.size, sizeof(*req_send));  // NOLINT

    irecv(j_recv, u, req_recv, nu);
    isend(i_send, send, u, buf, req_send, nu);
    wait(req_recv);
    wait(req_send);
}

void sync_begin(Equations *eqns, double *u, long nu)
{
    const ALIAS(j_recv, eqns->mesh->sync.j_recv);
    irecv(j_recv, u, eqns->sync.recv, nu);

    const ALIAS(i_send, eqns->mesh->sync.i_send);
    const ALIAS(send, eqns->mesh->sync.send);
    ALIAS(buf, eqns->sync.buf);
    isend(i_send, send, u, buf, eqns->sync.send, nu);
}

void sync_wait(Equations *eqns)
{
    wait(eqns->sync.recv);
}

void sync_end(Equations *eqns)
{
    wait(eqns->sync.send);
}

static void irecv(const long *j_recv, double *u, MPI_Request *req, long nu)
{
    for (long rank = 0; rank < teal.size; ++rank) {
        const long count = (j_recv[rank + 1] - j_recv[rank]) * nu;
        MPI_Irecv(&u[j_recv[rank] * nu], count, MPI_DOUBLE, rank, 0, teal.comm, &req[rank]);
    }
}

static void isend(const long *i_send, const long *send, const double *u, double *buf,
                  MPI_Request *req, long nu)
{
    for (long rank = 0; rank < teal.size; ++rank)
        for (long i = i_send[rank]; i < i_send[rank + 1]; ++i)
            for (long j = 0; j < nu; ++j) buf[i * nu + j] = u[send[i] * nu + j];

    for (long rank = 0; rank < teal.size; ++rank) {
        const long count = (i_send[rank + 1] - i_send[rank]) * nu;
        MPI_Isend(&buf[i_send[rank] * nu], count, MPI_DOUBLE, rank, 0, teal.comm, &req[rank]);
    }
}

static void wait(MPI_Request *req)
{
    MPI_Waitall(teal.size, req, MPI_STATUSES_IGNORE);
}
