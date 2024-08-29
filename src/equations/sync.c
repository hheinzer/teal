#include "sync.h"

#include "teal/sync.h"

void sync_begin(Equations *eqns, double *u, long ldu)
{
    sync_irecv(eqns->mesh->sync.j_recv, eqns->sync.recv, u, ldu);
    sync_isend(eqns->mesh->sync.i_send, eqns->mesh->sync.send, u, eqns->sync.send, eqns->sync.buf,
               ldu);
}

void sync_wait(Equations *eqns)
{
    sync_waitall(eqns->sync.recv);
}

void sync_end(Equations *eqns)
{
    sync_waitall(eqns->sync.send);
}
