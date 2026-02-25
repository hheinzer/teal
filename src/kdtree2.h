#pragma once

#include "vector2.h"

typedef struct kdtree Kdtree2;

// Initialize a KD-tree that stores fixed-size payloads of `size` bytes.
Kdtree2 *kdtree2_init(int size);

// Free a KD-tree and all stored payload copies.
void kdtree2_deinit(Kdtree2 *self);

// Insert a point and copy `size` bytes from `data` into the KD-tree.
void kdtree2_insert(Kdtree2 *self, Vector pos, const void *data);

// Copy up to `cap` nearest payloads into `data`; returns number copied.
int kdtree2_nearest(const Kdtree2 *self, Vector pos, void *data, int cap);

// Copy up to `cap` payloads within `radius`; returns total matches (may exceed `cap`).
int kdtree2_radius(const Kdtree2 *self, Vector pos, void *data, int cap, double radius);
