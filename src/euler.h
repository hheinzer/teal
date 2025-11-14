#pragma once

#include "equations.h"
#include "simulation.h"  // IWYU pragma: export

typedef struct {
    scalar density;
    vector momentum;
    scalar energy;
    vector velocity;
    scalar pressure;
} Euler;

Update euler_conserved, euler_primitive;
Timestep euler_timestep;
BoundarySelect euler_boundary;
ConvectiveSelect euler_convective;

Equations *euler_create(const Mesh *mesh);

/* Solve the exact Riemann problem for the Euler equations. */
Euler euler_riemann(const Euler *left, const Euler *right, scalar gamma, scalar location);
