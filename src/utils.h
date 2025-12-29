#pragma once

#define sq(val) ((val) * (val))
#define cb(val) ((val) * (val) * (val))

#define min(lhs, rhs) (((lhs) < (rhs)) ? (lhs) : (rhs))
#define max(lhs, rhs) (((lhs) > (rhs)) ? (lhs) : (rhs))

int isclose(double lhs, double rhs);

int cmp_long(const void *lhs_, const void *rhs_);
int cmp_double(const void *lhs_, const void *rhs_);
