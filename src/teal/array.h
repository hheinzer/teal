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
vector array_vwsum(const vector *arr, const double *wgt, long num);

long array_ldigitize(const long *arr, long val, long num);
long array_fdigitize(const double *arr, double val, long num);
