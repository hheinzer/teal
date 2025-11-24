#pragma once

#include "teal.h"

long array_lmin(const long *arr, long num);
long array_lmax(const long *arr, long num);
long array_lsum(const long *arr, long num);

scalar array_fmin(const scalar *arr, long num);
scalar array_fmax(const scalar *arr, long num);
scalar array_fsum(const scalar *arr, long num);

void array_unique(long *arr, long *num);

/* Return index of first `arr[i] > val` in sorted array; if beyond bounds return `0` or `num`. */
long array_digitize(const long *arr, long val, long num);

scalar array_dot(const scalar *lhs, const scalar *rhs, long num);
scalar array_norm(const scalar *arr, long num);
