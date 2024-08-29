#pragma once

#include "teal.h"

void euler_global_to_local(const Matrix3d b, const double *g, double *l);

void euler_local_to_global(const Matrix3d b, const double *l, double *g);
