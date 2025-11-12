#include "sync.h"

#include "teal/assert.h"
#include "teal/sync.h"

Request sync_variables(const Equations *eqns, void *variable, number stride)
{
    assert(eqns && variable && stride > 0);

    number num = eqns->mesh->neighbors.num;
    number *rank = eqns->mesh->neighbors.rank;
    number *off_recv = eqns->mesh->neighbors.recv_off;
    number *off_send = eqns->mesh->neighbors.send.off;
    number *idx_send = eqns->mesh->neighbors.send.idx;

    int tag = sync_tag();

    Request req;
    req.recv = sync_irecv_scalar(rank, off_recv, variable, num, stride, tag);
    req.send = sync_isend_scalar(rank, off_send, idx_send, variable, num, stride, tag);

    return req;
}

Request sync_gradients(const Equations *eqns, void *gradient, number stride)
{
    assert(eqns && gradient && stride > 0);

    number num = eqns->mesh->neighbors.num;
    number *rank = eqns->mesh->neighbors.rank;
    number *off_recv = eqns->mesh->neighbors.recv_off;
    number *off_send = eqns->mesh->neighbors.send.off;
    number *idx_send = eqns->mesh->neighbors.send.idx;

    int tag = sync_tag();

    Request req;
    req.recv = sync_irecv_vector(rank, off_recv, gradient, num, stride, tag);
    req.send = sync_isend_vector(rank, off_send, idx_send, gradient, num, stride, tag);

    return req;
}

void sync_wait(const Equations *eqns, MPI_Request *req)
{
    assert(eqns && req);
    number num = eqns->mesh->neighbors.num;
    MPI_Waitall(num, req, MPI_STATUSES_IGNORE);
}
