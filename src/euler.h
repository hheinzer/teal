#ifndef EULER_H
#define EULER_H

#include "equations.h"
#include "mesh.h"

#ifndef VARIABLES
#define VARIABLES
typedef enum Variables {
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
} Variables;
#endif

#ifndef SCALARS
#define SCALARS
typedef enum Scalars {
    GAMMA,  // heat capacity ratio
    N_SCALARS,
} Scalars;
#endif

Equations euler_create(const Mesh *mesh, long space_order);

#endif
