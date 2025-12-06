#pragma once

#include "equations.h"
#include "simulation.h"  // IWYU pragma: export

typedef struct {
    scalar density;
    vector momentum;
    scalar energy;
    vector velocity;
    scalar pressure;
} NavierStokes;

typedef enum {
    NAVIER_STOKES_HEAT_CAPACITY_RATIO,
    NAVIER_STOKES_DYNAMIC_VISCOSITY,
    NAVIER_STOKES_PRANDTL,
} NavierStokesProperties;

// Compute the Navier-Stokes time step restriction.
TimeStep navier_stokes_time_step;

// Select a Navier-Stokes boundary condition by name.
BoundarySelect navier_stokes_boundary;

// Select a viscous Navier-Stokes flux implementation.
ViscousSelect navier_stokes_viscous;

// Create the Navier-Stokes equations on a mesh with default properties and fluxes.
Equations *navier_stokes_create(const Mesh *mesh);
