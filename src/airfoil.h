#ifndef AIRFOIL_H
#define AIRFOIL_H

#include "mesh.h"
#include "simulation.h"

Mesh airfoil_mesh(const char *fname, const long n_upper, const long n_lower, const long n_outer,
                  const double radius);

void airfoil_polar(const Simulation *sim, const char *name, const long P, const double *uref,
                   const double lref);

#endif
