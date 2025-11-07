#include "array.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "assert.h"
#include "utils.h"

number array_lmin(const number *arr, number num)
{
    assert(arr ? num >= 0 : num == 0);
    number min = PTRDIFF_MAX;
    for (number i = 0; i < num; i++) {
        min = lmin(min, arr[i]);
    }
    return min;
}

number array_lmax(const number *arr, number num)
{
    assert(arr ? num >= 0 : num == 0);
    number max = PTRDIFF_MIN;
    for (number i = 0; i < num; i++) {
        max = lmax(max, arr[i]);
    }
    return max;
}

number array_lsum(const number *arr, number num)
{
    assert(arr ? num >= 0 : num == 0);
    number sum = 0;
    for (number i = 0; i < num; i++) {
        sum += arr[i];
    }
    return sum;
}

scalar array_fmin(const scalar *arr, number num)
{
    assert(arr ? num >= 0 : num == 0);
    scalar min = INFINITY;
    for (number i = 0; i < num; i++) {
        min = fmin(min, arr[i]);
    }
    return min;
}

scalar array_fmax(const scalar *arr, number num)
{
    assert(arr ? num >= 0 : num == 0);
    scalar max = -INFINITY;
    for (number i = 0; i < num; i++) {
        max = fmax(max, arr[i]);
    }
    return max;
}

scalar array_fsum(const scalar *arr, number num)
{
    assert(arr ? num >= 0 : num == 0);
    double sum = 0;  // always use double to avoid loss of precision
    for (number i = 0; i < num; i++) {
        sum += arr[i];
    }
    return sum;
}

void array_lunique(number *arr, number *num)
{
    assert(arr ? *num >= 0 : *num == 0);
    if (*num <= 1) {
        return;
    }
    qsort(arr, *num, sizeof(*arr), cmp_number);
    number unique = 1;
    for (number i = 1; i < *num; i++) {
        if (arr[i] != arr[unique - 1]) {
            arr[unique++] = arr[i];
        }
    }
    *num = unique;
}

number array_ldigitize(const number *arr, number val, number num)
{
    assert(arr ? num >= 0 : num == 0);
    if (num == 0) {
        return 0;
    }
    number left = 0;
    number right = num - 1;
    if (val < arr[left]) {
        return 0;
    }
    if (arr[right] <= val) {
        return num;
    }
    while (left < right) {
        number mid = left + ((right - left) / 2);
        if (arr[mid] <= val) {
            left = mid + 1;
        }
        else {
            right = mid;
        }
    }
    return left;
}
