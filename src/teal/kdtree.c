#include "kdtree.h"

#include <assert.h>
#include <float.h>
#include <stdio.h>

#include "isclose.h"
#include "memory.h"

static void insert(Kdtree *tree, const double *key, const long *val, long nval, KdtreeItem *root,
                   long depth);

static void nearest(const Kdtree *tree, const double *key, KdtreeItem *item, KdtreeItem *root,
                    long depth, long count);

static void radius(const Kdtree *tree, const double *key, KdtreeItem **item, KdtreeItem *root,
                   long depth, long *count, double r2);

static double metric(const double *a, const double *b, long n);

static void replace(KdtreeItem *heap, KdtreeItem *new, double dist2, long count);

static void append(KdtreeItem **item, KdtreeItem *new, double dist2, long *count);

Kdtree *kdtree_create(long max_items, long nkey)
{
    Kdtree *tree = memory_calloc(1, sizeof(*tree));
    tree->n_items = 0;
    tree->max_items = max_items;
    tree->nkey = nkey;
    tree->item = memory_calloc(tree->max_items, sizeof(*tree->item));
    return tree;
}

void kdtree_insert(Kdtree *tree, const double *key, const long *val, long nval)
{
    assert(tree->n_items < tree->max_items);
    if (tree->n_items == 0) {
        tree->item->key = key;
        tree->item->nval = nval;
        tree->item->val = memory_duplicate(val, nval, sizeof(*val));
        tree->n_items += 1;
    }
    else
        insert(tree, key, val, nval, tree->item, 0);
}

KdtreeItem *kdtree_nearest(const Kdtree *tree, const double *key, long count)
{
    assert(count <= tree->n_items);
    KdtreeItem *item = memory_calloc(count, sizeof(*item));
    for (long i = 0; i < count; ++i) item[i].dist2 = DBL_MAX;
    nearest(tree, key, item, tree->item, 0, count);
    return item;
}

KdtreeItem *kdtree_radius(const Kdtree *tree, const double *key, long *count, double r)
{
    *count = 0;
    KdtreeItem *item = 0;
    radius(tree, key, &item, tree->item, 0, count, r * r);
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

void kdtree_free(Kdtree **tree)
{
    for (long i = 0; i < (*tree)->n_items; ++i) memory_free(&(*tree)->item[i].val);
    memory_free(&(*tree)->item);
    memory_free(tree);
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
        assert(!is_close(metric(key, root->key, tree->nkey), 0));
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
                    long depth, long count)
{
    if (!root) return;
    const long axis = depth % tree->nkey;
    const double dx = key[axis] - root->key[axis];
    const double dist2 = metric(key, root->key, tree->nkey);
    if (dist2 < item->dist2) replace(item, root, dist2, count);
    nearest(tree, key, item, (dx < 0 ? root->left : root->right), depth + 1, count);
    if (dx * dx < item->dist2)
        nearest(tree, key, item, (dx < 0 ? root->right : root->left), depth + 1, count);
}

static void radius(const Kdtree *tree, const double *key, KdtreeItem **item, KdtreeItem *root,
                   long depth, long *count, double r2)
{
    if (!root) return;
    const long axis = depth % tree->nkey;
    const double dx = key[axis] - root->key[axis];
    const double dist2 = metric(key, root->key, tree->nkey);
    if (dist2 <= r2) append(item, root, dist2, count);
    radius(tree, key, item, (dx < 0 ? root->left : root->right), depth + 1, count, r2);
    if (dx * dx < r2)
        radius(tree, key, item, (dx < 0 ? root->right : root->left), depth + 1, count, r2);
}

static double metric(const double *a, const double *b, long n)
{
    double m = 0;
    for (long i = 0; i < n; ++i) m += (a[i] - b[i]) * (a[i] - b[i]);
    return m;
}

static void replace(KdtreeItem *heap, KdtreeItem *new, double dist2, long count)
{
    heap->key = new->key;
    heap->nval = new->nval;
    heap->val = new->val;
    heap->dist2 = dist2;
    long i = 0, left = 1;
    while (left < count) {
        const long right = left + 1;
        long max = i;
        if (heap[left].dist2 > heap[max].dist2) max = left;
        if (right < count && heap[right].dist2 > heap[max].dist2) max = right;
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

static void append(KdtreeItem **item, KdtreeItem *new, double dist2, long *count)
{
    *item = memory_realloc(*item, *count + 1, sizeof(**item));
    (*item)[*count].key = new->key;
    (*item)[*count].nval = new->nval;
    (*item)[*count].val = new->val;
    (*item)[*count].dist2 = dist2;
    *count += 1;
}
