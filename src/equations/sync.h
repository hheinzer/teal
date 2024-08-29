#pragma once

#include "equations.h"

void sync_begin(Equations *eqns, double *u, long ldu);

void sync_wait(Equations *eqns);

void sync_end(Equations *eqns);
