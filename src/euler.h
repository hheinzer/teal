#ifndef EULER_H
#define EULER_H

#include "equations.h"
#include "mesh.h"

#ifndef VARIABLES
#define VARIABLES
enum Variables {
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

#ifndef SCALARS
#define SCALARS
enum Scalars {
    GAMMA,  // heat capacity ratio
    N_SCALARS,
};
#endif

Equations euler_create(const Mesh *mesh, long space_order);

#endif
