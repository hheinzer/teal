#include "sync.h"

#include <assert.h>
#include <mpi.h>

#include "equations.h"
#include "memory.h"
#include "mesh.h"

static void irecv(const int size, const long *i_recv, const long nu, double (*u)[nu],
                  MPI_Request *req);
static void isend(const int size, const long *i_send, const long *send, const long nu,
                  const double (*u)[nu], double (*buf)[nu], MPI_Request *req);
static void wait(const int size, MPI_Request *req);

long sync_max_long(const long v)
{
    long r = 0;
    MPI_Allreduce(&v, &r, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
    return r;
}
long *sync_max_long_n(long *v, const long n)
{
    MPI_Allreduce(MPI_IN_PLACE, v, n, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
    return v;
}
double sync_max_double(const double v)
{
    double r = 0;
    MPI_Allreduce(&v, &r, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return r;
}
double *sync_max_double_n(double *v, const long n)
{
    MPI_Allreduce(MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return v;
}

long sync_min_long(const long v)
{
    long r = 0;
    MPI_Allreduce(&v, &r, 1, MPI_LONG, MPI_MIN, MPI_COMM_WORLD);
    return r;
}
long *sync_min_long_n(long *v, const long n)
{
    MPI_Allreduce(MPI_IN_PLACE, v, n, MPI_LONG, MPI_MIN, MPI_COMM_WORLD);
    return v;
}
double sync_min_double(const double v)
{
    double r = 0;
    MPI_Allreduce(&v, &r, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return r;
}
double *sync_min_double_n(double *v, const long n)
{
    MPI_Allreduce(MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return v;
}

long sync_sum_long(const long v)
{
    long r = 0;
    MPI_Allreduce(&v, &r, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    return r;
}
long *sync_sum_long_n(long *v, const long n)
{
    MPI_Allreduce(MPI_IN_PLACE, v, n, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    return v;
}
double sync_sum_double(const double v)
{
    double r = 0;
    MPI_Allreduce(&v, &r, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return r;
}
double *sync_sum_double_n(double *v, const long n)
{
    MPI_Allreduce(MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return v;
}

long sync_exscan_sum_long(const long v)
{
    long r = 0;
    MPI_Exscan(&v, &r, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    return r;
}
long *sync_exscan_sum_long_n(long *v, const long n)
{
    MPI_Exscan(MPI_IN_PLACE, v, n, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    return v;
}
double sync_exscan_sum_double(const double v)
{
    double r = 0;
    MPI_Exscan(&v, &r, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return r;
}
double *sync_exscan_sum_double_n(double *v, const long n)
{
    MPI_Exscan(MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return v;
}

void sync_all(const Mesh *mesh, const long nu, double (*u)[nu])
{
    const int size = mesh->size;
    const ALIAS(i_recv, mesh->sync.i_recv);
    const ALIAS(i_send, mesh->sync.i_send);
    const ALIAS(send, mesh->sync.send);

    cleanup double(*buf)[nu] = memory_calloc(i_send[size], sizeof(*buf));
    cleanup MPI_Request(*req)[size] = memory_calloc(2, sizeof(*req));

    irecv(size, i_recv, nu, u, req[0]);
    isend(size, i_send, send, nu, u, buf, req[1]);
    wait(2 * size, *req);
}

void sync_u_begin(Equations *eqns)
{
    const int size = eqns->mesh->size;
    const ALIAS(i_recv, eqns->mesh->sync.i_recv);
    const long nu = eqns->vars.n_fields;
    FIELDS(u, eqns->vars);
    irecv(size, i_recv, nu, u, eqns->sync.recv_u);

    const ALIAS(i_send, eqns->mesh->sync.i_send);
    const ALIAS(send, eqns->mesh->sync.send);
    double(*buf)[nu] = (typeof(buf))eqns->sync.buf_u;
    isend(size, i_send, send, nu, u, buf, eqns->sync.send_u);
}

void sync_u_wait(Equations *eqns) { wait(eqns->mesh->size, eqns->sync.recv_u); }

void sync_u_end(Equations *eqns) { wait(eqns->mesh->size, eqns->sync.send_u); }

void sync_dudx_begin(Equations *eqns)
{
    const int size = eqns->mesh->size;
    const ALIAS(i_recv, eqns->mesh->sync.i_recv);
    const long nu = eqns->vars.n_fields * N_DIMS;
    GRADS(dudx, eqns->vars);
    irecv(size, i_recv, nu, *dudx, eqns->sync.recv_dudx);

    const ALIAS(i_send, eqns->mesh->sync.i_send);
    const ALIAS(send, eqns->mesh->sync.send);
    double(*buf)[nu] = (typeof(buf))eqns->sync.buf_dudx;
    isend(size, i_send, send, nu, *dudx, buf, eqns->sync.send_dudx);
}

void sync_dudx_wait(Equations *eqns) { wait(eqns->mesh->size, eqns->sync.recv_dudx); }

void sync_dudx_end(Equations *eqns) { wait(eqns->mesh->size, eqns->sync.send_dudx); }

static void irecv(const int size, const long *i_recv, const long nu, double (*u)[nu],
                  MPI_Request *req)
{
    for (long rank = 0; rank < size; ++rank) {
        const long count = (i_recv[rank + 1] - i_recv[rank]) * nu;
        const int status =
            MPI_Irecv(&u[i_recv[rank]], count, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD, &req[rank]);
        assert(status == MPI_SUCCESS && "irecv error");
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
        const int status =
            MPI_Isend(buf[i_send[rank]], count, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD, &req[rank]);
        assert(status == MPI_SUCCESS && "isend error");
    }
}

static void wait(const int size, MPI_Request *req)
{
    const int status = MPI_Waitall(size, req, MPI_STATUSES_IGNORE);
    assert(status == MPI_SUCCESS && "wait error");
}
