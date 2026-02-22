#pragma once

#include <parmetis.h>

#include "grid.h"

_Static_assert(sizeof(idx_t) == sizeof(long), "idx_t must match long");
_Static_assert(sizeof(real_t) == sizeof(double), "real_t must match double");

typedef struct {
    long *dist;
    long *xadj;
    long *adjncy;
    long *part;
} Dual;

Dual *dual_init(const Grid *grid);

void dual_deinit(Dual *dual);

void dual_partition(Dual *dual, const Grid *grid);

void dual_periodic(Dual *dual, const Grid *grid);

int compare_idx(const void *lhs, const void *rhs);
