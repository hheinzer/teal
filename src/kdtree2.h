#pragma once

#include "vector2.h"

typedef struct kdtree Kdtree2;

// Build a KD-tree over `num` points.
Kdtree2 *kdtree2_init(const Vector *point, int num);

// Free a KD-tree created with `tree2_init`.
void kdtree2_deinit(Kdtree2 *self);

// Copy up to `cap` nearest point indices into `idx`; returns number copied.
int kdtree2_nearest(const Kdtree2 *self, Vector point, int *idx, int cap);

// Copy up to `cap` point indices within `radius`; returns total matches.
int kdtree2_radius(const Kdtree2 *self, Vector point, double radius, int *idx, int cap);
