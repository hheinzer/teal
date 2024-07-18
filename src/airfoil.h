#ifndef AIRFOIL_H
#define AIRFOIL_H

#include "simulation.h"

void airfoil_polar(const Simulation *sim, const char *airfoil, const double *uref, double lref,
                   long D, long U, long V, long W, long P);

#endif
