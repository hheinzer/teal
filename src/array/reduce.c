#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include "array.h"
#include "utils.h"

int array_min(const int *arr, int num)
{
    assert((arr || num == 0) && num >= 0);
    int min_val = INT_MAX;
    for (int i = 0; i < num; i++) {
        min_val = min(min_val, arr[i]);
    }
    return min_val;
}

int array_max(const int *arr, int num)
{
    assert((arr || num == 0) && num >= 0);
    int max_val = INT_MIN;
    for (int i = 0; i < num; i++) {
        max_val = max(max_val, arr[i]);
    }
    return max_val;
}

int array_sum(const int *arr, int num)
{
    assert((arr || num == 0) && num >= 0);
    int sum = 0;
    for (int i = 0; i < num; i++) {
        sum += arr[i];
    }
    return sum;
}

double array_fmin(const double *arr, int num)
{
    assert((arr || num == 0) && num >= 0);
    double min_val = DBL_MAX;
    for (int i = 0; i < num; i++) {
        min_val = fmin(min_val, arr[i]);
    }
    return min_val;
}

double array_fmax(const double *arr, int num)
{
    assert((arr || num == 0) && num >= 0);
    double max_val = -DBL_MAX;
    for (int i = 0; i < num; i++) {
        max_val = fmax(max_val, arr[i]);
    }
    return max_val;
}

double array_fsum(const double *arr, int num)
{
    assert((arr || num == 0) && num >= 0);
    double sum = 0;
    for (int i = 0; i < num; i++) {
        sum += arr[i];
    }
    return sum;
}
