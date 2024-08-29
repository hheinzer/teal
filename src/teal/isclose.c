#include "isclose.h"

#include <math.h>

// https://numpy.org/doc/stable/reference/generated/numpy.isclose.html#numpy-isclose
#define ATOL 1e-8
#define RTOL 1e-5

int is_close(double a, double b)
{
    // https://docs.python.org/3/library/math.html#math.isclose
    const double tol = fmax(ATOL, RTOL * fmax(fabs(a), fabs(b)));
    return fabs(a - b) <= tol;
}

int is_close_or_less(double a, double b)
{
    return is_close(a, b) || a < b;
}

int is_close_or_greater(double a, double b)
{
    return is_close(a, b) || a > b;
}

int is_less(double a, double b)
{
    return !is_close(a, b) && a < b;
}

int is_greater(double a, double b)
{
    return !is_close(a, b) && a > b;
}
