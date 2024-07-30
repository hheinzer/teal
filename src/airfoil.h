#pragma once

#include "simulation.h"

/* Compute airfoil polar of entity 'airfoil' using reference velocity 'uref' and reference length
 * 'lref'. 'D', 'U', 'V', 'W', and 'P' are the indices of density, x-, y-, z-velocity, and pressure
 * of the equation system in 'sim'. */
void airfoil_polar(const Simulation *sim, const char *airfoil, const double *uref, double lref,
                   long D, long U, long V, long W, long P);
