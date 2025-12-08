#include "sync.h"

#include <assert.h>

#include "teal/sync.h"

Request sync_variables(const Equations *eqns, void *variable, long stride)
{
    assert(eqns && variable && stride > 0);

    long num = eqns->mesh->neighbors.num;
    long *rank = eqns->mesh->neighbors.rank;
    long *off_recv = eqns->mesh->neighbors.recv_off;
    long *off_send = eqns->mesh->neighbors.send.off;
    long *idx_send = eqns->mesh->neighbors.send.idx;

    long tag = sync_tag();

    Request req;
    req.recv = sync_irecv(rank, off_recv, variable, num, stride, MPI_SCALAR, tag);
    req.send = sync_isend(rank, off_send, idx_send, variable, num, stride, MPI_SCALAR, tag);
    return req;
}

Request sync_gradients(const Equations *eqns, void *gradient, long stride)
{
    assert(eqns && gradient && stride > 0);

    long num = eqns->mesh->neighbors.num;
    long *rank = eqns->mesh->neighbors.rank;
    long *off_recv = eqns->mesh->neighbors.recv_off;
    long *off_send = eqns->mesh->neighbors.send.off;
    long *idx_send = eqns->mesh->neighbors.send.idx;

    long tag = sync_tag();

    Request req;
    req.recv = sync_irecv(rank, off_recv, gradient, num, stride * 3, MPI_SCALAR, tag);
    req.send = sync_isend(rank, off_send, idx_send, gradient, num, stride * 3, MPI_SCALAR, tag);
    return req;
}

void sync_wait(const Equations *eqns, MPI_Request *req)
{
    assert(eqns && req);
    long num = eqns->mesh->neighbors.num;
    double wtime_beg = MPI_Wtime();
    MPI_Waitall(num, req, MPI_STATUSES_IGNORE);
    double wtime_end = MPI_Wtime();
    sync.wait += wtime_end - wtime_beg;
}
