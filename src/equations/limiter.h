#pragma once

#include "equations.h"

Limiter vanleer, vanalbada1, vanalbada2, mc, koren, minmod, superbee, venkatakrishnan;

// Precompute per-cell epsilon^2 for the Venkatakrishnan limiter using cell volumes.
scalar *venkatakrishnan_parameter(const Equations *eqns, scalar parameter);
