#include <assert.h>
#include <stdlib.h>

#include "array.h"
#include "utils.h"

void array_unique(long *arr, long *num)
{
    assert((arr || *num == 0) && *num >= 0);
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

long array_digitize(const long *arr, long val, long num)
{
    assert((arr || num == 0) && num >= 0);
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
