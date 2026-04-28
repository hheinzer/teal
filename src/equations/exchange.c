#include "exchange.h"

#include <string.h>

#include "sync.h"
#include "teal.h"

Exchange equations_exchange(const Equations *eqns, void *buf_, int len)
{
    int num = eqns->mesh->neighbors.num;
    int (*tag)[2] = eqns->mesh->neighbors.tag;
    int *rank = eqns->mesh->neighbors.rank;
    int *recv_off = eqns->mesh->neighbors.recv_off;
    int *send_off = eqns->mesh->neighbors.send.off;
    int *send_idx = eqns->mesh->neighbors.send.idx;

    double (*buf)[len] = buf_;

    double (*send)[len] = teal_calloc(send_off[num], (int)sizeof(*send));
    for (int i = 0; i < num; i++) {
        for (int j = send_off[i]; j < send_off[i + 1]; j++) {
            memcpy(send[j], buf[send_idx[j]], sizeof(*buf));
        }
    }

    MPI_Request *req_recv = teal_calloc(num, sizeof(*req_recv));
    MPI_Request *req_send = teal_calloc(num, sizeof(*req_send));
    for (int i = 0; i < num; i++) {
        int num_recv = recv_off[i + 1] - recv_off[i];
        int num_send = send_off[i + 1] - send_off[i];
        sync_irecv(&req_recv[i], buf[recv_off[i]], num_recv, rank[i], tag[i][0], MPI_DOUBLE, len);
        sync_isend(&req_send[i], send[send_off[i]], num_send, rank[i], tag[i][1], MPI_DOUBLE, len);
    }

    return (Exchange){send, req_recv, req_send};
}

void equations_exchange_wait_recv(const Equations *eqns, Exchange exchange)
{
    int num = eqns->mesh->neighbors.num;
    MPI_Waitall(num, exchange.req_recv, MPI_STATUSES_IGNORE);
    teal_free(exchange.req_recv);
}

void equations_exchange_wait_send(const Equations *eqns, Exchange exchange)
{
    int num = eqns->mesh->neighbors.num;
    MPI_Waitall(num, exchange.req_send, MPI_STATUSES_IGNORE);
    teal_free(exchange.req_send);
    teal_free(exchange.send);
}
