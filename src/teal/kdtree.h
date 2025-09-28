#pragma once

#include "teal.h"

typedef struct Kdtree Kdtree;
typedef struct KdtreeItem KdtreeItem;

struct Kdtree {
    long num;
    KdtreeItem *beg;
    KdtreeItem **end;
    long size_val;
};

struct KdtreeItem {
    vector key;
    void *val;
    KdtreeItem *left;
    KdtreeItem *right;
    KdtreeItem *next;
};

/* Returns empty kd-tree; if `size_val == 0` it behaves as a set of points (values ignored). */
Kdtree *kdtree_create(long size_val);

/* Returns pointer to stored value if `key` is already present; else inserts and returns `0`. */
void *kdtree_insert(Kdtree *self, vector key, const void *val);

/* Returns pointer to stored value if `key` is present or `0` if absent. */
void *kdtree_lookup(const Kdtree *self, vector key);

/* Writes up to `num` nearest values to `key` into `val` (unsorted). */
void kdtree_nearest(const Kdtree *self, vector key, void *val, long num);

/* Writes up to `cap` values within `radius` to `key` into `val` (unsorted); returns count. */
long kdtree_radius(const Kdtree *self, vector key, void *val, long cap, double radius);
