#pragma once

#include "euler.h"

/* Solve the exact Riemann problem for the Euler equations. */
Euler riemann(const Euler *left, const Euler *right, scalar gamma, scalar location);
