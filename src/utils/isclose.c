#include <math.h>

#include "utils.h"

static const double tol_absolute = 1e-8;
static const double tol_relative = 1e-5;

int isclose(double lhs, double rhs)
{
    if (lhs == rhs) {
        return 1;
    }
    if (!isfinite(lhs) || !isfinite(rhs)) {
        return 0;
    }
    double diff = fabs(lhs - rhs);
    if (diff <= tol_absolute) {
        return 1;
    }
    return diff <= tol_relative * fmax(fabs(lhs), fabs(rhs));
}
