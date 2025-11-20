#include "array.h"

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"

int array_lmin(const int *arr, int num)
{
    assert(arr ? (num >= 0) : (num == 0));
    int min = INT_MAX;
    for (int i = 0; i < num; i++) {
        min = lmin(min, arr[i]);
    }
    return min;
}

int array_lmax(const int *arr, int num)
{
    assert(arr ? (num >= 0) : (num == 0));
    int max = INT_MIN;
    for (int i = 0; i < num; i++) {
        max = lmax(max, arr[i]);
    }
    return max;
}

int array_lsum(const int *arr, int num)
{
    assert(arr ? (num >= 0) : (num == 0));
    int sum = 0;
    for (int i = 0; i < num; i++) {
        sum += arr[i];
    }
    return sum;
}

scalar array_fmin(const scalar *arr, int num)
{
    assert(arr ? (num >= 0) : (num == 0));
    scalar min = INFINITY;
    for (int i = 0; i < num; i++) {
        min = fmin(min, arr[i]);
    }
    return min;
}

scalar array_fmax(const scalar *arr, int num)
{
    assert(arr ? (num >= 0) : (num == 0));
    scalar max = -INFINITY;
    for (int i = 0; i < num; i++) {
        max = fmax(max, arr[i]);
    }
    return max;
}

scalar array_fsum(const scalar *arr, int num)
{
    assert(arr ? (num >= 0) : (num == 0));
    double sum = 0;  // always use double to avoid loss of precision
    for (int i = 0; i < num; i++) {
        sum += arr[i];
    }
    return sum;
}

void array_lunique(int *arr, int *num)
{
    assert(arr ? (*num >= 0) : (*num == 0));
    if (*num <= 1) {
        return;
    }
    qsort(arr, *num, sizeof(*arr), cmp_int);
    int unique = 1;
    for (int i = 1; i < *num; i++) {
        if (arr[i] != arr[unique - 1]) {
            arr[unique++] = arr[i];
        }
    }
    *num = unique;
}

int array_ldigitize(const int *arr, int val, int num)
{
    assert(arr ? (num >= 0) : (num == 0));
    if (num == 0) {
        return 0;
    }
    int left = 0;
    int right = num - 1;
    if (val < arr[left]) {
        return 0;
    }
    if (arr[right] <= val) {
        return num;
    }
    while (left < right) {
        int mid = left + ((right - left) / 2);
        if (arr[mid] <= val) {
            left = mid + 1;
        }
        else {
            right = mid;
        }
    }
    return left;
}

scalar array_fdot(const scalar *lhs, const scalar *rhs, int num)
{
    assert((lhs || rhs) ? (num >= 0) : (num == 0));
    double dot = 0;
    for (int i = 0; i < num; i++) {
        dot += lhs[i] * rhs[i];
    }
    return dot;
}

scalar array_fnorm(const scalar *arr, int num)
{
    return sqrt(array_fdot(arr, arr, num));
}
