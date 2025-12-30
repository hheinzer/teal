#include <math.h>

#include "array.h"
#include "sync.h"

double sync_dot(const double *lhs, const double *rhs, int num)
{
    return sync_fsum(array_dot(lhs, rhs, num));
}

double sync_norm(const double *arr, int num)
{
    return sqrt(sync_dot(arr, arr, num));
}
