#pragma once

long array_lmin(const long *arr, long num);
long array_lmax(const long *arr, long num);
long array_lsum(const long *arr, long num);

double array_fmin(const double *arr, long num);
double array_fmax(const double *arr, long num);
double array_fsum(const double *arr, long num);

// Sort the array and remove duplicates in place.
void array_unique(long *arr, long *num);

// Return index of first element greater than val in a sorted array.
long array_digitize(const long *arr, long val, long num);

// Compute the dot product of two arrays.
double array_dot(const double *lhs, const double *rhs, long num);

// Compute the Euclidean norm of an array.
double array_norm(const double *arr, long num);
