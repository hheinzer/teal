#include "array.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>

#include "assert.h"
#include "utils.h"
#include "vector.h"

long array_lmin(const long *arr, long num)
{
    assert(arr ? num >= 0 : num == 0);
    long min = LONG_MAX;
    for (long i = 0; i < num; i++) {
        min = lmin(min, arr[i]);
    }
    return min;
}

long array_lmax(const long *arr, long num)
{
    assert(arr ? num >= 0 : num == 0);
    long max = LONG_MIN;
    for (long i = 0; i < num; i++) {
        max = lmax(max, arr[i]);
    }
    return max;
}

long array_lsum(const long *arr, long num)
{
    assert(arr ? num >= 0 : num == 0);
    long sum = 0;
    for (long i = 0; i < num; i++) {
        sum += arr[i];
    }
    return sum;
}

scalar array_fmin(const scalar *arr, long num)
{
    assert(arr ? num >= 0 : num == 0);
    scalar min = INFINITY;
    for (long i = 0; i < num; i++) {
        min = fmin(min, arr[i]);
    }
    return min;
}

scalar array_fmax(const scalar *arr, long num)
{
    assert(arr ? num >= 0 : num == 0);
    scalar max = -INFINITY;
    for (long i = 0; i < num; i++) {
        max = fmax(max, arr[i]);
    }
    return max;
}

scalar array_fsum(const scalar *arr, long num)
{
    assert(arr ? num >= 0 : num == 0);
    double sum = 0;  // always use double to avoid loss of precision
    for (long i = 0; i < num; i++) {
        sum += arr[i];
    }
    return sum;
}

vector array_vmin(const vector *arr, long num)
{
    assert(arr ? num >= 0 : num == 0);
    vector min = {INFINITY, INFINITY, INFINITY};
    for (long i = 0; i < num; i++) {
        min = vector_min(min, arr[i]);
    }
    return min;
}

vector array_vmax(const vector *arr, long num)
{
    assert(arr ? num >= 0 : num == 0);
    vector max = {-INFINITY, -INFINITY, -INFINITY};
    for (long i = 0; i < num; i++) {
        max = vector_max(max, arr[i]);
    }
    return max;
}

vector array_vsum(const vector *arr, long num)
{
    assert(arr ? num >= 0 : num == 0);
    vector sum = {0};
    for (long i = 0; i < num; i++) {
        sum = vector_add(sum, arr[i]);
    }
    return sum;
}

void array_lunique(long *arr, long *num)
{
    assert(arr ? *num >= 0 : *num == 0);
    if (*num <= 1) {
        return;
    }
    qsort(arr, *num, sizeof(*arr), cmp_long);
    long unique = 1;
    for (long i = 1; i < *num; i++) {
        if (arr[i] != arr[unique - 1]) {
            arr[unique++] = arr[i];
        }
    }
    *num = unique;
}

long array_ldigitize(const long *arr, long val, long num)
{
    assert(arr ? num >= 0 : num == 0);
    if (num == 0) {
        return 0;
    }
    long left = 0;
    long right = num - 1;
    if (val < arr[left]) {
        return 0;
    }
    if (arr[right] <= val) {
        return num;
    }
    while (left < right) {
        long mid = left + ((right - left) / 2);
        if (arr[mid] <= val) {
            left = mid + 1;
        }
        else {
            right = mid;
        }
    }
    return left;
}
