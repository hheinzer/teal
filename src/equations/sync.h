#pragma once

#include <mpi.h>

#include "equations.h"

typedef struct {
    MPI_Request *recv;
    MPI_Request *send;
} Request;

Request sync_variables(const Equations *eqns, void *variable, int stride);

Request sync_gradients(const Equations *eqns, void *gradient, int stride);

void sync_wait(const Equations *eqns, MPI_Request *req);
