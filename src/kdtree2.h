#pragma once

#include "vector2.h"

typedef struct kdtree Kdtree;

// Build a KD-tree over `num` points.
Kdtree *kdtree2_init(const Vector *point, int num);

// Free a KD-tree created with `tree2_init`.
void kdtree2_deinit(Kdtree *self);

// Copy up to `cap` nearest point indices into `idx`; returns number copied.
int kdtree2_nearest(const Kdtree *self, Vector query, int *idx, int cap);

// Copy up to `cap` point indices within `radius`; returns total matches.
int kdtree2_radius(const Kdtree *self, Vector query, double radius, int *idx, int cap);
