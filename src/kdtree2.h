#pragma once

#include "vector2.h"

typedef struct kdtree Kdtree;

// Build a KD-tree over `num` points.
Kdtree *kdtree2_init(const Vector *point, int num);

// Free all memory associated with the tree.
void kdtree2_deinit(Kdtree *self);

// Copy up to `cap` nearest point indices into `index`; returns number copied.
int kdtree2_nearest(const Kdtree *self, Vector point, int *index, int cap);
