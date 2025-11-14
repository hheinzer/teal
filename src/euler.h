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

BoundarySelect euler_boundary;
ConvectiveSelect euler_convective;

Equations *euler_create(const Mesh *mesh);
