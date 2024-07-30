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

/* Create a k-d tree with enough space of 'max_items' and a key length of 'nkey'. */
Kdtree kdtree_create(long max_items, long nkey);

/* Deallocate memory of 'tree' and set all fields to 0. */
void kdtree_free(Kdtree *tree);

/* Insert 'key'-'val' pair with value length 'nval' into 'tree'. Collisions are not allowed. */
void kdtree_insert(Kdtree *tree, const double *key, const long *val, long nval);

/* Find the 'k' nearest items to 'key' in 'tree'. */
KdtreeItem *kdtree_nearest(const Kdtree *tree, const double *key, long k);

/* Find the 'k' items that have a radius of 'r' or less to 'key' in 'tree'. */
KdtreeItem *kdtree_radius(const Kdtree *tree, const double *key, long *k, double r);

/* Print content of 'tree'. */
void kdtree_print(const Kdtree *tree);
