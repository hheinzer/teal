#pragma once

#include <mpi.h>

#include "equations.h"

typedef struct {
    void *send;
    MPI_Request *req_recv;
    MPI_Request *req_send;
} Exchange;

Exchange equations_exchange(const Equations *eqns, void *buf, int len);

void equations_exchange_wait_recv(const Equations *eqns, Exchange exchange);

void equations_exchange_wait_send(const Equations *eqns, Exchange exchange);
