#pragma once

#include <stddef.h>

#define sq(val) ((val) * (val))
#define cb(val) ((val) * (val) * (val))

#define min(lhs, rhs) (((lhs) < (rhs)) ? (lhs) : (rhs))
#define max(lhs, rhs) (((lhs) > (rhs)) ? (lhs) : (rhs))

#define inrange(min, val, max) (min <= val && val <= max)

// Return nonzero if two doubles are close within tolerances.
int isclose(double lhs, double rhs);

// Compare two ints for qsort ordering.
int cmp_int(const void *lhs_, const void *rhs_);

// Compare two longs for qsort ordering.
int cmp_long(const void *lhs_, const void *rhs_);

// Compare two doubles for qsort ordering.
int cmp_double(const void *lhs_, const void *rhs_);

void *memdup(const void *ptr, int num, size_t size);
