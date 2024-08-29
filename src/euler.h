#pragma once

#include "equations.h"
#include "mesh.h"
#include "simulation.h"

enum Conserved { D, DU, DV, DW, DE, N_CONS };
enum Primitive { U = N_CONS, V, W, P, N_VARS };
enum Scalar { GAMMA, N_SCALARS };

Equations *euler_create(const Mesh *mesh);

void euler_riemann(double *d, double *u, double *p, double dl, double ul, double pl, double dr,
                   double ur, double pr, double s, double gamma);

void euler_polar(const Simulation *sim, const char *airfoil, const double *uref, double lref);
