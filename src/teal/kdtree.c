#include "kdtree.h"

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "arena.h"
#include "utils.h"
#include "vector.h"

#define NON_NULL ((void *)(int)1)  // NOLINT(performance-no-int-to-ptr)

Kdtree *kdtree_create(int size_val)
{
    assert(size_val >= 0);
    Kdtree *tree = arena_calloc(1, sizeof(*tree));
    tree->end = &tree->beg;
    tree->size_val = size_val;
    return tree;
}

static int veccmp(vector lhs, vector rhs, int depth)
{
    switch (depth % 3) {
        case 0: return isclose(lhs.x, rhs.x) ? 0 : cmp_asc(lhs.x, rhs.x);
        case 1: return isclose(lhs.y, rhs.y) ? 0 : cmp_asc(lhs.y, rhs.y);
        case 2: return isclose(lhs.z, rhs.z) ? 0 : cmp_asc(lhs.z, rhs.z);
        default: abort();
    }
}

static int keycmp(vector lhs, vector rhs, int depth)
{
    int cmp = veccmp(lhs, rhs, depth);
    if (cmp) {
        return cmp;
    }
    cmp = veccmp(lhs, rhs, depth + 1);
    if (cmp) {
        return cmp;
    }
    return veccmp(lhs, rhs, depth + 2);
}

void *kdtree_insert(Kdtree *self, vector key, const void *val)
{
    assert(self && (val || self->size_val == 0));

    KdtreeItem **item = &self->beg;
    int depth = 0;
    while (*item) {
        int cmp = keycmp(key, (*item)->key, depth);
        if (!cmp) {
            return (*item)->val;
        }
        item = (cmp < 0) ? &(*item)->left : &(*item)->right;
        depth += 1;
    }

    *item = arena_calloc(1, sizeof(**item));
    (*item)->key = key;
    (*item)->val = self->size_val ? arena_memdup(val, 1, self->size_val) : NON_NULL;

    self->num += 1;
    *self->end = *item;
    self->end = &(*item)->next;

    return 0;
}

void *kdtree_lookup(const Kdtree *self, vector key)
{
    assert(self);
    KdtreeItem *item = self->beg;
    int depth = 0;
    while (item) {
        int cmp = keycmp(key, item->key, depth);
        if (!cmp) {
            return item->val;
        }
        item = (cmp < 0) ? item->left : item->right;
        depth += 1;
    }
    return 0;
}

static scalar squared_distance(vector lhs, vector rhs)
{
    vector sub = vector_sub(lhs, rhs);
    return vector_dot(sub, sub);
}

static void heap_replace(const void *item_val, scalar item_metric, void *val_, scalar *metric,
                         int num, int size)
{
    char (*val)[size] = val_;

    memcpy(val[0], item_val, size);
    metric[0] = item_metric;

    int idx = 0;
    int left = 1;
    while (left < num) {
        int right = left + 1;
        int max = idx;
        if (metric[left] > metric[max]) {
            max = left;
        }
        if (right < num && metric[right] > metric[max]) {
            max = right;
        }
        if (max != idx) {
            swap_bytes(val[idx], val[max], size);
            swap_scalar(&metric[idx], &metric[max]);
            idx = max;
            left = (2 * idx) + 1;
        }
        else {
            break;
        }
    }
}

static scalar delta(vector lhs, vector rhs, int depth)
{
    switch (depth % 3) {
        case 0: return lhs.x - rhs.x;
        case 1: return lhs.y - rhs.y;
        case 2: return lhs.z - rhs.z;
        default: abort();
    }
}

typedef struct {
    KdtreeItem *item;
    int depth;
} Stack;

void kdtree_nearest(const Kdtree *self, vector key, void *val, int num)
{
    assert(self && self->beg && self->size_val > 0 && val && 0 < num && num <= self->num);
    Arena save = arena_save();

    scalar *metric = arena_malloc(num, sizeof(*metric));
    for (int i = 0; i < num; i++) {
        metric[i] = INFINITY;
    }

    Stack *stack = arena_malloc(self->num, sizeof(*stack));
    int top = 0;

    stack[top++] = (Stack){self->beg, 0};
    while (top) {
        Stack cur = stack[--top];

        scalar item_metric = squared_distance(key, cur.item->key);
        if (item_metric < metric[0]) {
            heap_replace(cur.item->val, item_metric, val, metric, num, self->size_val);
        }

        scalar del = delta(key, cur.item->key, cur.depth);
        KdtreeItem *near = (del < 0) ? cur.item->left : cur.item->right;
        KdtreeItem *far = (del < 0) ? cur.item->right : cur.item->left;

        // push far first so near is popped next (improve pruning)
        if (far && sq(del) <= metric[0]) {
            stack[top++] = (Stack){far, cur.depth + 1};
        }
        if (near) {
            stack[top++] = (Stack){near, cur.depth + 1};
        }
    }

    arena_load(save);
}

int kdtree_radius(const Kdtree *self, vector key, void *val_, int cap, scalar radius)
{
    assert(self && self->beg && self->size_val > 0 && val_ && 0 < cap && radius >= 0);
    Arena save = arena_save();

    char (*val)[self->size_val] = val_;
    scalar metric = sq(radius);
    int num = 0;

    Stack *stack = arena_malloc(self->num, sizeof(*stack));
    int top = 0;

    stack[top++] = (Stack){self->beg, 0};
    while (top && num < cap) {
        Stack cur = stack[--top];

        scalar item_metric = squared_distance(key, cur.item->key);
        if (item_metric <= metric) {
            memcpy(val[num++], cur.item->val, self->size_val);
        }

        scalar del = delta(key, cur.item->key, cur.depth);
        KdtreeItem *near = (del < 0) ? cur.item->left : cur.item->right;
        KdtreeItem *far = (del < 0) ? cur.item->right : cur.item->left;

        // push far first so near is popped next (improve pruning)
        if (far && sq(del) <= metric) {
            stack[top++] = (Stack){far, cur.depth + 1};
        }
        if (near) {
            stack[top++] = (Stack){near, cur.depth + 1};
        }
    }

    arena_load(save);
    return num;
}
