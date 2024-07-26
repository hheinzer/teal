#pragma once

#include "equations.h"
#include "mesh.h"

#ifndef VARIABLE
#define VARIABLE
enum Variable {
    D,   // density
    U,   // velocity-x
    V,   // velocity-y
    W,   // velocity-z
    P,   // pressure
    DU,  // momentum-x
    DV,  // momentum-y
    DW,  // momentum-z
    DE,  // energy
    N_VARS,
};
#endif

#ifndef SCALAR
#define SCALAR
enum Scalar {
    GAMMA,    // heat capacity ratio
    PRANDTL,  // Prandtl number
    MU,       // dynamic viscosity
    N_SCALARS,
};
#endif

Equations navierstokes_create(const Mesh *mesh, long space_order);
