#pragma once

#include "equations.h"
#include "simulation.h"  // IWYU pragma: export
#include "teal.h"        // IWYU pragma: export

typedef struct {
    double density;
    Vector velocity;
    double pressure;
} NavierStokesPrimitive;

typedef struct {
    double density;
    Vector momentum;
    double energy;
} NavierStokesConserved;

enum NavierStokesProperty {
    NAVIERSTOKES_HEAT_CAPACITY_RATIO,
    NAVIERSTOKES_DYNAMIC_VISCOSITY,
    NAVIERSTOKES_PRANDTL,
};

Timestep navierstokes_timestep;
ViscousSelect navierstokes_viscous;
BoundarySelect navierstokes_boundary;

// Create an Navier-Stokes equation system.
Equations *navierstokes_create(const Mesh *mesh);
