#pragma once

#include "equations.h"
#include "simulation.h"
#include "teal.h"  // IWYU pragma: export

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

Timestep euler_timestep;
ConvectiveSelect euler_convective;
BoundarySelect euler_boundary;
Convert euler_primitive, euler_conserved;

// Create an Euler equation system.
Equations *euler_create(const Mesh *mesh);

// Solve the exact Riemann problem for the Euler equations.
EulerPrimitive euler_riemann(EulerPrimitive left, EulerPrimitive right, double gamma, double loc);

// Compute pressure coefficient distribution and lift/drag coefficients for an airfoil entity.
void euler_polar(const Simulation *sim, const char *entity, EulerPrimitive reference, double length,
                 double time);
