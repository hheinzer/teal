#pragma once

#include "equations.h"
#include "mesh.h"
#include "simulation.h"

enum Conserved { D, DU, DV, DW, DE, N_CONS };
enum Primitive { U = N_CONS, V, W, P, N_VARS };
enum Scalar { GAMMA, MU, PRANDTL, N_SCALARS };

Equations *navierstokes_create(const Mesh *mesh);
