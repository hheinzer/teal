#include "sync.h"

#include "teal/assert.h"
#include "teal/sync.h"

Request sync_variables(const Equations *eqns, void *variable, int stride)
{
    assert(eqns && variable && stride > 0);

    int num = eqns->mesh->neighbors.num;
    int *rank = eqns->mesh->neighbors.rank;
    int *off_recv = eqns->mesh->neighbors.recv_off;
    int *off_send = eqns->mesh->neighbors.send.off;
    int *idx_send = eqns->mesh->neighbors.send.idx;

    int tag = sync_tag();

    Request req;
    req.recv = sync_irecv(rank, off_recv, variable, num, stride, MPI_SCALAR, tag);
    req.send = sync_isend(rank, off_send, idx_send, variable, num, stride, MPI_SCALAR, tag);
    return req;
}

Request sync_gradients(const Equations *eqns, void *gradient, int stride)
{
    assert(eqns && gradient && stride > 0);

    int num = eqns->mesh->neighbors.num;
    int *rank = eqns->mesh->neighbors.rank;
    int *off_recv = eqns->mesh->neighbors.recv_off;
    int *off_send = eqns->mesh->neighbors.send.off;
    int *idx_send = eqns->mesh->neighbors.send.idx;

    int tag = sync_tag();

    Request req;
    req.recv = sync_irecv(rank, off_recv, gradient, num, stride * 3, MPI_SCALAR, tag);
    req.send = sync_isend(rank, off_send, idx_send, gradient, num, stride * 3, MPI_SCALAR, tag);
    return req;
}

void sync_wait(const Equations *eqns, MPI_Request *req)
{
    assert(eqns && req);
    int num = eqns->mesh->neighbors.num;
    MPI_Waitall(num, req, MPI_STATUSES_IGNORE);
}
