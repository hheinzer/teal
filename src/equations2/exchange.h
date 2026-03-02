#pragma once

#include <mpi.h>

#include "equations2.h"

typedef struct {
    void *send;
    MPI_Request *req_recv;
    MPI_Request *req_send;
} Exchange;

Exchange equations2_exchange(const Equations *eqns, void *buf, int len);

void equations2_exchange_wait_recv(const Equations *eqns, Exchange exchange);

void equations2_exchange_wait_send(const Equations *eqns, Exchange exchange);
