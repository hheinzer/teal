#pragma once

#include "equations2.h"

typedef struct {
    double density;
    Vector velocity;
    double pressure;
} EulerPrimitive;

typedef struct {
    double density;
    Vector momentum;
    double energy;
} EulerConserved;

enum EulerProperty {
    EULER_HEAT_CAPACITY_RATIO,
};

Timestep euler2_timestep;
ConvectiveSelect euler2_convective;
BoundarySelect euler2_boundary;
Convert euler2_primitive, euler2_conserved;

// Create an Euler equation system.
Equations *euler2_create(const Mesh *mesh);

// Solve the exact Riemann problem for the Euler equations.
EulerPrimitive euler2_riemann(EulerPrimitive left, EulerPrimitive right, double gamma, double loc);
