#pragma once

#include "teal.h"

long array_lmin(const long *arr, long num);
long array_lmax(const long *arr, long num);
long array_lsum(const long *arr, long num);

double array_fmin(const double *arr, long num);
double array_fmax(const double *arr, long num);
double array_fsum(const double *arr, long num);

vector array_vmin(const vector *arr, long num);
vector array_vmax(const vector *arr, long num);
vector array_vsum(const vector *arr, long num);

/* Return index of first `arr[i] > val` in sorted array; if beyond bounds return `0` or `num`. */
long array_ldigitize(const long *arr, long val, long num);
