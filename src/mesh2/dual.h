#pragma once

#include <parmetis.h>

#include "grid.h"

typedef struct {
    idx_t *dist;
    idx_t *xadj;
    idx_t *adjncy;
    idx_t *part;
} Dual;

Dual *dual_init(const Grid *grid);

void dual_deinit(Dual *dual);

void dual_partition(Dual *dual, const Grid *grid);

void dual_periodic(Dual *dual, const Grid *grid);

int compare_idx(const void *lhs, const void *rhs);
