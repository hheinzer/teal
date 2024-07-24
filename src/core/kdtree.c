#include "kdtree.h"

#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#include "core/memory.h"
#include "core/utils.h"

static void insert(Kdtree *tree, const double *key, const long *val, long nval, KdtreeItem *root,
                   long depth);

static void nearest(const Kdtree *tree, const double *key, KdtreeItem *item, KdtreeItem *root,
                    long depth, long k);

static void radius(const Kdtree *tree, const double *key, KdtreeItem **item, KdtreeItem *root,
                   long depth, long *k, double r2);

static double metric(const double *a, const double *b, long n);

static void replace(KdtreeItem *heap, KdtreeItem *new, double dist2, long size);

static void append(KdtreeItem **item, KdtreeItem *new, double dist2, long *size);

Kdtree kdtree_create(long max_items, long nkey)
{
    Kdtree tree = {
        .n_items = 0,
        .max_items = max_items,
        .nkey = nkey,
    };
    tree.item = memory_calloc(tree.max_items, sizeof(*tree.item));
    return tree;
}

void kdtree_free(Kdtree *tree)
{
    for (long i = 0; i < tree->n_items; ++i) free(tree->item[i].val);
    free(tree->item);
    *tree = (Kdtree){0};
}

void kdtree_insert(Kdtree *tree, const double *key, const long *val, long nval)
{
    ensure(tree->n_items < tree->max_items);
    if (tree->n_items == 0) {
        tree->item->key = key;
        tree->item->nval = nval;
        tree->item->val = memory_duplicate(val, nval, sizeof(*val));
        tree->n_items += 1;
    }
    else
        insert(tree, key, val, nval, tree->item, 0);
}

KdtreeItem *kdtree_nearest(const Kdtree *tree, const double *key, long k)
{
    ensure(k <= tree->n_items);
    KdtreeItem *item = memory_calloc(k, sizeof(*item));
    for (long i = 0; i < k; ++i) item[i].dist2 = DBL_MAX;
    nearest(tree, key, item, tree->item, 0, k);
    return item;
}

KdtreeItem *kdtree_radius(const Kdtree *tree, const double *key, long *k, double r)
{
    *k = 0;
    KdtreeItem *item = 0;
    radius(tree, key, &item, tree->item, 0, k, sq(r));
    return item;
}

void kdtree_print(const Kdtree *tree)
{
    for (long i = 0; i < tree->n_items; ++i) {
        printf("%ld: ", i);
        for (long j = 0; j < tree->nkey; ++j) printf("%g ", tree->item[i].key[j]);
        printf(": ");
        for (long j = 0; j < tree->item[i].nval; ++j) printf("%ld ", tree->item[i].val[j]);
        printf(": ");
        printf("%td ", (tree->item[i].left ? tree->item[i].left - &tree->item[0] : 0));
        printf("%td ", (tree->item[i].right ? tree->item[i].right - &tree->item[0] : 0));
        printf("\n");
    }
}

static void insert(Kdtree *tree, const double *key, const long *val, long nval, KdtreeItem *root,
                   long depth)
{
    const long axis = depth % tree->nkey;
    const double dx = key[axis] - root->key[axis];
    KdtreeItem *next = (dx < 0 ? root->left : root->right);
    if (next)
        insert(tree, key, val, nval, next, depth + 1);
    else {
        if (dx < 0)
            next = root->left = &tree->item[tree->n_items];
        else
            next = root->right = &tree->item[tree->n_items];
        next->key = key;
        next->nval = nval;
        next->val = memory_duplicate(val, nval, sizeof(*val));
        tree->n_items += 1;
    }
}

static void nearest(const Kdtree *tree, const double *key, KdtreeItem *item, KdtreeItem *root,
                    long depth, long k)
{
    if (!root) return;
    const long axis = depth % tree->nkey;
    const double dx = key[axis] - root->key[axis];
    const double dist2 = metric(key, root->key, tree->nkey);
    if (dist2 < item->dist2) replace(item, root, dist2, k);
    nearest(tree, key, item, (dx < 0 ? root->left : root->right), depth + 1, k);
    if (sq(dx) < item->dist2)
        nearest(tree, key, item, (dx < 0 ? root->right : root->left), depth + 1, k);
}

static void radius(const Kdtree *tree, const double *key, KdtreeItem **item, KdtreeItem *root,
                   long depth, long *k, double r2)
{
    if (!root) return;
    const long axis = depth % tree->nkey;
    const double dx = key[axis] - root->key[axis];
    const double dist2 = metric(key, root->key, tree->nkey);
    if (dist2 <= r2) append(item, root, dist2, k);
    radius(tree, key, item, (dx < 0 ? root->left : root->right), depth + 1, k, r2);
    if (sq(dx) < r2) radius(tree, key, item, (dx < 0 ? root->right : root->left), depth + 1, k, r2);
}

static double metric(const double *a, const double *b, long n)
{
    double m = 0;
    for (long i = 0; i < n; ++i) m += sq(a[i] - b[i]);
    return m;
}

static void replace(KdtreeItem *heap, KdtreeItem *new, double dist2, long size)
{
    heap->key = new->key;
    heap->nval = new->nval;
    heap->val = new->val;
    heap->dist2 = dist2;
    long i = 0, left = 1;
    while (left < size) {
        const long right = left + 1;
        long max = i;
        if (heap[left].dist2 > heap[max].dist2) max = left;
        if (right < size && heap[right].dist2 > heap[max].dist2) max = right;
        if (max != i) {
            KdtreeItem swap = heap[i];
            heap[i] = heap[max];
            heap[max] = swap;
            i = max;
            left = 2 * i + 1;
        }
        else
            break;
    }
}

static void append(KdtreeItem **item, KdtreeItem *new, double dist2, long *size)
{
    *item = memory_realloc(*item, *size + 1, sizeof(**item));
    (*item)[*size].key = new->key;
    (*item)[*size].nval = new->nval;
    (*item)[*size].val = new->val;
    (*item)[*size].dist2 = dist2;
    *size += 1;
}
