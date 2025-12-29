#include <math.h>

#include "array.h"
#include "sync.h"

double sync_dot(const double *lhs, const double *rhs, long num)
{
    return sync_fsum(array_dot(lhs, rhs, num));
}

double sync_norm(const double *arr, long num)
{
    return sqrt(sync_dot(arr, arr, num));
}
