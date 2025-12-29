#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include "array.h"
#include "utils.h"

long array_lmin(const long *arr, long num)
{
    assert((arr || num == 0) && num >= 0);
    long min_val = LONG_MAX;
    for (long i = 0; i < num; i++) {
        min_val = min(min_val, arr[i]);
    }
    return min_val;
}

long array_lmax(const long *arr, long num)
{
    assert((arr || num == 0) && num >= 0);
    long max_val = LONG_MIN;
    for (long i = 0; i < num; i++) {
        max_val = max(max_val, arr[i]);
    }
    return max_val;
}

long array_lsum(const long *arr, long num)
{
    assert((arr || num == 0) && num >= 0);
    long sum = 0;
    for (long i = 0; i < num; i++) {
        sum += arr[i];
    }
    return sum;
}

double array_fmin(const double *arr, long num)
{
    assert((arr || num == 0) && num >= 0);
    double min_val = DBL_MAX;
    for (long i = 0; i < num; i++) {
        min_val = fmin(min_val, arr[i]);
    }
    return min_val;
}

double array_fmax(const double *arr, long num)
{
    assert((arr || num == 0) && num >= 0);
    double max_val = -DBL_MAX;
    for (long i = 0; i < num; i++) {
        max_val = fmax(max_val, arr[i]);
    }
    return max_val;
}

double array_fsum(const double *arr, long num)
{
    assert((arr || num == 0) && num >= 0);
    double sum = 0;
    for (long i = 0; i < num; i++) {
        sum += arr[i];
    }
    return sum;
}
