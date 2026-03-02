#pragma once

#include "equations2.h"

enum LimiterKind {
    NONE,
    BARTH_JESPERSEN,
    VENKATAKRISHNAN,
};

double *venkatakrishnan_parameter2(const Equations *eqns, double coefficient);
