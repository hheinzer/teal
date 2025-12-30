#include <assert.h>
#include <float.h>
#include <math.h>

#include "vector.h"

vector vector_min(const vector *arr, int num)
{
    assert((arr || num == 0) && num >= 0);
    vector min = {.x = DBL_MAX, .y = DBL_MAX, .z = DBL_MAX};
    for (int i = 0; i < num; i++) {
        min.x = fmin(min.x, arr[i].x);
        min.y = fmin(min.y, arr[i].y);
        min.z = fmin(min.z, arr[i].z);
    }
    return min;
}

vector vector_max(const vector *arr, int num)
{
    assert((arr || num == 0) && num >= 0);
    vector max = {.x = -DBL_MAX, .y = -DBL_MAX, .z = -DBL_MAX};
    for (int i = 0; i < num; i++) {
        max.x = fmax(max.x, arr[i].x);
        max.y = fmax(max.y, arr[i].y);
        max.z = fmax(max.z, arr[i].z);
    }
    return max;
}

vector vector_sum(const vector *arr, int num)
{
    assert((arr || num == 0) && num >= 0);
    vector sum = {0};
    for (int i = 0; i < num; i++) {
        sum.x += arr[i].x;
        sum.y += arr[i].y;
        sum.z += arr[i].z;
    }
    return sum;
}
