#pragma once

#include "simulation.h"

void airfoil_polar(const Simulation *sim, const char *airfoil, const double *uref, double lref,
                   long D, long U, long V, long W, long P);
