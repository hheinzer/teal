#include <assert.h>

#include "utils.h"

int cmp_int(const void *lhs_, const void *rhs_)
{
    assert(lhs_ && rhs_);
    const int *lhs = lhs_;
    const int *rhs = rhs_;
    return (*lhs > *rhs) - (*lhs < *rhs);
}

int cmp_double(const void *lhs_, const void *rhs_)
{
    assert(lhs_ && rhs_);
    const double *lhs = lhs_;
    const double *rhs = rhs_;
    return (*lhs > *rhs) - (*lhs < *rhs);
}
