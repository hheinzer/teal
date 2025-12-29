#include <assert.h>
#include <math.h>

#include "array.h"

double array_dot(const double *lhs, const double *rhs, long num)
{
    assert(((lhs && rhs) || num == 0) && num >= 0);
    double dot = 0;
    for (long i = 0; i < num; i++) {
        dot += lhs[i] * rhs[i];
    }
    return dot;
}

double array_norm(const double *arr, long num)
{
    assert((arr || num == 0) && num >= 0);
    return sqrt(array_dot(arr, arr, num));
}
