#include <assert.h>
#include <stdlib.h>

#include "array.h"
#include "utils.h"

void array_unique(int *arr, int *num)
{
    assert((arr || *num == 0) && *num >= 0);
    if (*num <= 1) {
        return;
    }
    qsort(arr, (size_t)*num, sizeof(*arr), cmp_int);
    int unique = 1;
    for (int i = 1; i < *num; i++) {
        if (arr[i] != arr[unique - 1]) {
            arr[unique++] = arr[i];
        }
    }
    *num = unique;
}

int array_digitize(const int *arr, int val, int num)
{
    assert((arr || num == 0) && num >= 0);
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
