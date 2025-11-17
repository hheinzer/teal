#pragma once

#include "simulation.h"

void gmres(const Equations *eqns, const void *variable_, const void *derivative_,
           const void *residual_, void *increment_, scalar time, scalar step, scalar norm,
           scalar fd_scale, const NewtonKrylov *ctx);
