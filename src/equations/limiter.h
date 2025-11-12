#pragma once

#include "equations.h"

Limiter minmod, venkatakrishnan;

scalar *venkatakrishnan_parameter(const Equations *eqns, scalar parameter);
