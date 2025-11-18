#pragma once

#include "teal.h"

int array_lmin(const int *arr, int num);
int array_lmax(const int *arr, int num);
int array_lsum(const int *arr, int num);

scalar array_fmin(const scalar *arr, int num);
scalar array_fmax(const scalar *arr, int num);
scalar array_fsum(const scalar *arr, int num);

void array_lunique(int *arr, int *num);

/* Return index of first `arr[i] > val` in sorted array; if beyond bounds return `0` or `num`. */
int array_ldigitize(const int *arr, int val, int num);

scalar array_fdot(const scalar *lhs, const scalar *rhs, int num);
scalar array_fnorm(const scalar *arr, int num);
