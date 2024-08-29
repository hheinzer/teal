#pragma once

#include "equations.h"

double navierstokes_timestep(const Equations *eqns, const double *u, const Vector3d projection,
                             double volume);
