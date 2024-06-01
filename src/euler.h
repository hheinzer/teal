#ifndef EULER_H
#define EULER_H

#include "equations.h"
#include "mesh.h"

typedef enum Variable {
    // primitive
    D,  // density (also conserved)
    U,  // x-velocity
    V,  // y-velocity
    P,  // pressure
    // conserved
    DU,  // x-momentum
    DV,  // y-momentum
    DE,  // energy
    N_VARS_EULER,
} Variable;

Equations euler_create(const Mesh *mesh, const Fields *user);

void euler_set_gamma(const double gamma);

void euler_compute_conserved(Equations *eqns);

void euler_compute_primitive(Equations *eqns);

void euler_print(const Equations *eqns);

Flux *euler_flux(const char *flux);

ApplyBC *euler_boundary_condition(const char *bc);

#endif
