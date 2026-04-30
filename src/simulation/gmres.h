#pragma once

#include "simulation.h"

void gmres(const Equations *eqns, const void *conserved, const void *derivative,
           const void *residual, void *increment, double time, double step, double norm,
           double fd_scale, const NewtonKrylov *ctx);
