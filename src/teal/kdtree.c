#include "kdtree.h"

#include <math.h>
#include <stddef.h>
#include <string.h>

#include "arena.h"
#include "utils.h"

Kdtree kdtree_create(long size_val)
{
    assert(size_val >= 0);
    Kdtree tree = {0};
    tree.size_val = size_val;
    return tree;
}

static int veccmp(vector lhs, vector rhs, long depth)
{
    switch (depth % 3) {
        case 0: return fcmp(&lhs.x, &rhs.x);
        case 1: return fcmp(&lhs.y, &rhs.y);
        case 2: return fcmp(&lhs.z, &rhs.z);
        default: assert(false);
    }
}

static int keycmp(vector lhs, vector rhs, long depth)
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
    long depth = 0;
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
    if (self->end) {
        self->end->next = *item;
    }
    self->end = *item;
    return 0;
}

void *kdtree_lookup(const Kdtree *self, vector key)
{
    assert(self);
    KdtreeItem *item = self->beg;
    long depth = 0;
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

typedef struct {
    KdtreeItem *item;
    long depth;
    void *prev;
} Stack;

static void push(Stack **self, KdtreeItem *item, long depth)
{
    Stack *next = arena_calloc(1, sizeof(*next));
    next->item = item;
    next->depth = depth;
    next->prev = *self;
    *self = next;
}

static Stack *pop(Stack **self)
{
    Stack *curr = *self;
    *self = curr->prev;
    return curr;
}

static double sqdist(vector lhs, vector rhs)
{
    vector sub = {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
    return (sub.x * sub.x) + (sub.y * sub.y) + (sub.z * sub.z);
}

static void heapreplace(const void *item_val, double item_metric, void *val, double *metric,
                        long num, long size)
{
    char (*_val)[size] = val;

    memcpy(_val[0], item_val, size);
    metric[0] = item_metric;

    long idx = 0;
    long left = 1;
    while (left < num) {
        long right = left + 1;
        long max = idx;
        if (metric[left] > metric[max]) {
            max = left;
        }
        if (right < num && metric[right] > metric[max]) {
            max = right;
        }
        if (max != idx) {
            memswap(_val[idx], _val[max], size);
            fswap(&metric[idx], &metric[max]);
            idx = max;
            left = 2 * idx + 1;
        }
        else {
            break;
        }
    }
}

static double delta(vector lhs, vector rhs, long depth)
{
    switch (depth % 3) {
        case 0: return lhs.x - rhs.x;
        case 1: return lhs.y - rhs.y;
        case 2: return lhs.z - rhs.z;
        default: assert(false);
    }
}

void kdtree_nearest(const Kdtree *self, vector key, void *val, long num)
{
    assert(self && (val ? (num > 0 && self->size_val > 0) : num == 0));
    if (num == 0) {
        return;
    }

    Arena save = arena_save();

    double *metric = arena_calloc(num, sizeof(*metric));
    for (long i = 0; i < num; i++) {
        metric[i] = INFINITY;
    }

    Stack *stack = 0;
    push(&stack, self->beg, 0);
    while (stack) {
        Stack *curr = pop(&stack);
        KdtreeItem *item = curr->item;
        long depth = curr->depth;
        if (!item) {
            continue;
        }

        double item_metric = sqdist(key, item->key);
        if (item_metric < metric[0]) {
            heapreplace(item->val, item_metric, val, metric, num, self->size_val);
        }

        double del = delta(key, item->key, depth);
        push(&stack, (del < 0) ? item->left : item->right, depth + 1);
        if (sq(del) <= metric[0]) {
            push(&stack, (del < 0) ? item->right : item->left, depth + 1);
        }
    }

    arena_load(save);
}

long kdtree_radius(const Kdtree *self, vector key, void *val, long cap, double radius)
{
    assert(self && (val ? (cap > 0 && self->size_val > 0) : cap == 0) && radius >= 0);
    if (cap == 0) {
        return 0;
    }

    Arena save = arena_save();

    char (*_val)[self->size_val] = val;
    double metric = sq(radius);
    long num = 0;

    Stack *stack = 0;
    push(&stack, self->beg, 0);
    while (stack && num < cap) {
        Stack *curr = pop(&stack);
        KdtreeItem *item = curr->item;
        long depth = curr->depth;
        if (!item) {
            continue;
        }

        double item_metric = sqdist(key, item->key);
        if (item_metric <= metric) {
            memcpy(_val[num++], item->val, self->size_val);
        }

        double del = delta(key, item->key, depth);
        push(&stack, (del < 0) ? item->left : item->right, depth + 1);
        if (sq(del) <= metric) {
            push(&stack, (del < 0) ? item->right : item->left, depth + 1);
        }
    }

    arena_load(save);
    return num;
}
