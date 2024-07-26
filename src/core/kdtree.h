#pragma once

typedef struct KdtreeItem {
    const double *key;
    long nval, *val;
    struct KdtreeItem *left, *right;
    double dist2;
} KdtreeItem;

typedef struct Kdtree {
    long n_items, max_items, nkey;
    KdtreeItem *item;
} Kdtree;

Kdtree kdtree_create(long max_items, long nkey);

void kdtree_free(Kdtree *tree);

void kdtree_insert(Kdtree *tree, const double *key, const long *val, long nval);

KdtreeItem *kdtree_nearest(const Kdtree *tree, const double *key, long k);

KdtreeItem *kdtree_radius(const Kdtree *tree, const double *key, long *k, double r);

void kdtree_print(const Kdtree *tree);
