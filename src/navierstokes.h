#pragma once

#include "equations.h"
#include "mesh.h"

#ifndef VARIABLE
#define VARIABLE
enum Variable { D, U, V, W, P, DU, DV, DW, DE, N_VARS };
#endif

#ifndef SCALAR
#define SCALAR
enum Scalar { GAMMA, PRANDTL, MU, N_SCALARS };
#endif

/* Create Navier-Stokes equation system on 'mesh' using a 'space_order' discretization. */
Equations navierstokes_create(const Mesh *mesh, long space_order);

/* Compute the body 'force' of 'eqns'. */
void navierstokes_body_force(const Equations *eqns, double *force);
