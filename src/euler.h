#pragma once

#include "equations.h"
#include "mesh.h"

#ifndef VARIABLE
#define VARIABLE
enum Conserved { D, DU, DV, DW, DE, N_CONS };
enum Primitive { U = N_CONS, V, W, P, N_VARS };
#endif

#ifndef SCALAR
#define SCALAR
enum Scalar { GAMMA, N_SCALARS };
#endif

/* Create Euler equation system on 'mesh' using a 'space_order' discretization. */
Equations euler_create(const Mesh *mesh, long space_order);
