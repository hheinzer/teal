#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "utils2.h"

int isclose(double lhs, double rhs)
{
    assert(isfinite(lhs) && isfinite(rhs));
    double abs_tol = 1e-8;
    double rel_tol = 1e-5;
    return fabs(lhs - rhs) <= fmax(abs_tol, rel_tol * fmax(fabs(lhs), fabs(rhs)));
}

int unique(void *base_, int num_, int size, Compare compare)
{
    assert((base_ || num_ == 0) && num_ >= 0 && size > 0 && compare);

    if (num_ <= 1) {
        return num_;
    }

    qsort(base_, num_, size, compare);

    char (*base)[size] = base_;

    int num = 1;
    for (int i = 1; i < num_; i++) {
        if (compare(base[i], base[num - 1]) != 0) {
            if (num != i) {
                memcpy(base[num], base[i], size);
            }
            num += 1;
        }
    }

    return num;
}

int digitize(const void *key, const void *base_, int num, int size, Compare compare)
{
    assert(key && (base_ || num == 0) && num >= 0 && size > 0 && compare);

    if (num == 0) {
        return -2;
    }

    const char (*base)[size] = base_;

    if (compare(key, base[0]) < 0) {
        return -1;
    }

    if (compare(key, base[1]) < 0) {
        return 0;
    }

    if (compare(key, base[num]) >= 0) {
        return num;
    }

    int left = 1;
    int right = num;
    while (left < right) {
        int middle = left + ((right - left) / 2);
        if (compare(base[middle], key) <= 0) {
            left = middle + 1;
        }
        else {
            right = middle;
        }
    }

    return left - 1;
}
