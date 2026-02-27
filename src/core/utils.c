#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "utils2.h"

int compare_int(const void *lhs_, const void *rhs_)
{
    const int *lhs = lhs_;
    const int *rhs = rhs_;
    return (*lhs > *rhs) - (*lhs < *rhs);
}

int compare_long(const void *lhs_, const void *rhs_)
{
    const long *lhs = lhs_;
    const long *rhs = rhs_;
    return (*lhs > *rhs) - (*lhs < *rhs);
}

int compare_double(const void *lhs_, const void *rhs_)
{
    const double *lhs = lhs_;
    const double *rhs = rhs_;
    return isgreater(*lhs, *rhs) - isless(*lhs, *rhs);
}

int isclose(double lhs, double rhs)
{
    static const double abs_tol = 1e-8;
    static const double rel_tol = 1e-5;
    return islessequal(fabs(lhs - rhs), fmax(abs_tol, rel_tol * fmax(fabs(lhs), fabs(rhs))));
}

void *copy(void *dst, const void *src, int num, int size)
{
    assert(((dst && src) || num == 0) && num >= 0 && size > 0);
    if (num == 0) {
        return 0;
    }
    return memcpy(dst, src, (size_t)num * size);
}

void *move(void *dst, const void *src, int num, int size)
{
    assert(((dst && src) || num == 0) && num >= 0 && size > 0);
    if (num == 0) {
        return 0;
    }
    return memmove(dst, src, (size_t)num * size);
}

void sort(void *base, int num, int size, Compare compare)
{
    assert((base || num == 0) && num >= 0 && size > 0 && compare);
    if (num <= 1) {
        return;
    }
    qsort(base, num, size, compare);
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

void *find(const void *key, const void *base_, int num, int size, Compare compare)
{
    assert(key && (base_ || num == 0) && num >= 0 && size > 0 && compare);

    if (num == 0) {
        return 0;
    }

    const char (*base)[size] = base_;

    for (int i = 0; i < num; i++) {
        if (compare(key, base[i]) == 0) {
            return (void *)base[i];
        }
    }

    return 0;
}

void *search(const void *key, const void *base, int num, int size, Compare compare)
{
    assert(key && (base || num == 0) && num >= 0 && size > 0 && compare);
    if (num == 0) {
        return 0;
    }
    return bsearch(key, base, num, size, compare);
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
