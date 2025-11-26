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

/* Convert primitive and conserved variables. */
Update euler_conserved, euler_primitive;

/* Compute the Euler time-step restriction. */
Timestep euler_timestep;

/* Select an Euler boundary condition by name. */
BoundarySelect euler_boundary;

/* Select a convective Euler flux implementation. */
ConvectiveSelect euler_convective;

/* Create the Euler equations on a mesh with default properties and fluxes. */
Equations *euler_create(const Mesh *mesh);

/* Solve the exact Riemann problem for the Euler equations. */
Euler euler_riemann(const Euler *left, const Euler *right, scalar gamma, scalar location);
