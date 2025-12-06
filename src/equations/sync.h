#pragma once

#include <mpi.h>

#include "equations.h"

typedef struct {
    MPI_Request *recv;
    MPI_Request *send;
} Request;

// Start nonblocking halo exchange for cell-centered variables.
Request sync_variables(const Equations *eqns, void *variable, long stride);

// Start nonblocking halo exchange for cell-centered gradients.
Request sync_gradients(const Equations *eqns, void *gradient, long stride);

// Wait for a previously initiated exchange.
void sync_wait(const Equations *eqns, MPI_Request *req);
