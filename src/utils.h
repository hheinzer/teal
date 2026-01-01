#pragma once

#include <stddef.h>

// Return the square of a value.
double sq(double val);

// Return the cube of a value.
double cb(double val);

// Return the minimum of two ints.
int min(int lhs, int rhs);

// Return the maximum of two ints.
int max(int lhs, int rhs);

// Return the minimum of two longs.
long lmin(long lhs, long rhs);

// Return the maximum of two longs.
long lmax(long lhs, long rhs);

// Return nonzero if two doubles are close within tolerances.
int isclose(double lhs, double rhs);

// Compare two ints for qsort ordering.
int cmp_int(const void *lhs_, const void *rhs_);

// Compare two longs for qsort ordering.
int cmp_long(const void *lhs_, const void *rhs_);

// Compare two doubles for qsort ordering.
int cmp_double(const void *lhs_, const void *rhs_);

// Duplicate num elements from a buffer.
void *memdup(const void *ptr, int num, size_t size);
