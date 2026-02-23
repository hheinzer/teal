#pragma once

#include "grid.h"

typedef struct {
    long *dist;
    long *xadj;
    long *adjncy;
} Dual;

Dual *dual_init(const Grid *grid);

void dual_deinit(Dual *dual);

long *dual_partition(Dual *dual, const Grid *grid);

void dual_periodic(Dual *dual, const Grid *grid);
