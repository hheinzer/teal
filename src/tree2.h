#pragma once

#include "vector2.h"

typedef struct tree Tree2;

// Build a KD-tree over `num` points with a leaf capacity limit.
Tree2 *tree2_init(const Vector *point, int num, int cap_leaf);

// Free a KD-tree created with `tree_init`.
void tree2_deinit(Tree2 *self);

// Copy up to `cap` nearest point indices into `idx`; returns number copied.
int tree2_nearest(const Tree2 *self, Vector point, int *idx, int cap);

// Copy up to `cap` point indices within `radius`; returns total matches.
int tree2_radius(const Tree2 *self, Vector point, double radius, int *idx, int cap);
