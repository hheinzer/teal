#pragma once

// Return the minimum of an integer array, or INT_MAX if num is 0.
int array_min(const int *arr, int num);

// Return the maximum of an integer array, or INT_MIN if num is 0.
int array_max(const int *arr, int num);

// Return the sum of an integer array.
int array_sum(const int *arr, int num);

// Return the minimum of a double array, or DBL_MAX if num is 0.
double array_fmin(const double *arr, int num);

// Return the maximum of a double array, or -DBL_MAX if num is 0.
double array_fmax(const double *arr, int num);

// Return the sum of a double array.
double array_fsum(const double *arr, int num);

// Sort the array and remove duplicates in place.
void array_unique(int *arr, int *num);

// Return index of first element greater than val in a sorted array.
int array_digitize(const int *arr, int val, int num);

// Compute the dot product of two arrays.
double array_dot(const double *lhs, const double *rhs, int num);

// Compute the Euclidean norm of an array.
double array_norm(const double *arr, int num);
