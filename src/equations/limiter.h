#pragma once

#include "equations.h"

enum LimiterKind {
    NONE,
    BARTH_JESPERSEN,
    VENKATAKRISHNAN,
};

double *venkatakrishnan_parameter(const Equations *eqns, double coefficient);
