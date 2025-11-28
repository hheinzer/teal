#pragma once

#include "equations.h"

// Limiter functions.
Limiter vanleer, vanalbada1, vanalbada2, mc, koren, minmod, superbee, venkatakrishnan;

// Precompute epsilon^2 for Venkatakrishnan limiter using cell volumes.
scalar *venkatakrishnan_parameter(const Equations *eqns, scalar parameter);
