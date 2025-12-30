#pragma once

#define sq(val) ((val) * (val))
#define cb(val) ((val) * (val) * (val))

#define min(lhs, rhs) (((lhs) < (rhs)) ? (lhs) : (rhs))
#define max(lhs, rhs) (((lhs) > (rhs)) ? (lhs) : (rhs))

#define in_range(min, val, max) (min <= val && val <= max)

// Return nonzero if two doubles are close within tolerances.
int isclose(double lhs, double rhs);

// Compare two ints for qsort ordering.
int cmp_int(const void *lhs_, const void *rhs_);

// Compare two doubles for qsort ordering.
int cmp_double(const void *lhs_, const void *rhs_);
