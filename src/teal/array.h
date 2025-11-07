#pragma once

#include "teal.h"

number array_lmin(const number *arr, number num);
number array_lmax(const number *arr, number num);
number array_lsum(const number *arr, number num);

scalar array_fmin(const scalar *arr, number num);
scalar array_fmax(const scalar *arr, number num);
scalar array_fsum(const scalar *arr, number num);

void array_lunique(number *arr, number *num);

/* Return index of first `arr[i] > val` in sorted array; if beyond bounds return `0` or `num`. */
number array_ldigitize(const number *arr, number val, number num);
