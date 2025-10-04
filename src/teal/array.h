#pragma once

#include "teal.h"

long array_lmin(const long *arr, long num);
long array_lmax(const long *arr, long num);
long array_lsum(const long *arr, long num);

scalar array_fmin(const scalar *arr, long num);
scalar array_fmax(const scalar *arr, long num);
scalar array_fsum(const scalar *arr, long num);

vector array_vmin(const vector *arr, long num);
vector array_vmax(const vector *arr, long num);
vector array_vsum(const vector *arr, long num);

void array_lunique(long *arr, long *num);

/* Return index of first `arr[i] > val` in sorted array; if beyond bounds return `0` or `num`. */
long array_ldigitize(const long *arr, long val, long num);
