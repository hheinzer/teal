#include "kdtree.h"

#include <math.h>
#include <stddef.h>
#include <string.h>

#include "arena.h"
#include "teal/vector.h"
#include "utils.h"

Kdtree *kdtree_create(long size_val)
{
    assert(size_val >= 0);
    Kdtree *tree = arena_calloc(1, sizeof(*tree));
    tree->end = &tree->beg;
    tree->size_val = size_val;
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
    if ((cmp = veccmp(lhs, rhs, depth + 1))) {
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
    *self->end = *item;
    self->end = &(*item)->next;

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

typedef struct Stack Stack;
struct Stack {
    KdtreeItem *item;
    long depth;
    Stack *prev;
};

static void push(Stack **top, KdtreeItem *item, long depth)
{
    Stack *new = arena_malloc(1, sizeof(*new));
    new->item = item;
    new->depth = depth;
    new->prev = *top;
    *top = new;
}

static Stack *pop(Stack **top)
{
    Stack *cur = *top;
    *top = cur->prev;
    return cur;
}

static scalar squared_distance(vector lhs, vector rhs)
{
    vector sub = vector_sub(lhs, rhs);
    return vector_dot(sub, sub);
}

static void heap_replace(const void *item_val, scalar item_metric, void *val, scalar *metric,
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

static scalar delta(vector lhs, vector rhs, long depth)
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

    scalar *metric = arena_malloc(num, sizeof(*metric));
    for (long i = 0; i < num; i++) {
        metric[i] = INFINITY;
    }

    Stack *top = 0;
    push(&top, self->beg, 0);
    while (top) {
        Stack *cur = pop(&top);
        KdtreeItem *item = cur->item;
        long depth = cur->depth;
        if (!item) {
            continue;
        }

        scalar item_metric = squared_distance(key, item->key);
        if (item_metric < metric[0]) {
            heap_replace(item->val, item_metric, val, metric, num, self->size_val);
        }

        scalar del = delta(key, item->key, depth);
        push(&top, (del < 0) ? item->left : item->right, depth + 1);
        if (sq(del) <= metric[0]) {
            push(&top, (del < 0) ? item->right : item->left, depth + 1);
        }
    }

    arena_load(save);
}

long kdtree_radius(const Kdtree *self, vector key, void *val, long cap, scalar radius)
{
    assert(self && (val ? (cap > 0 && self->size_val > 0) : cap == 0) && radius >= 0);
    if (cap == 0) {
        return 0;
    }

    Arena save = arena_save();

    char (*_val)[self->size_val] = val;
    scalar metric = sq(radius);
    long num = 0;

    Stack *top = 0;
    push(&top, self->beg, 0);
    while (top && num < cap) {
        Stack *cur = pop(&top);
        KdtreeItem *item = cur->item;
        long depth = cur->depth;
        if (!item) {
            continue;
        }

        scalar item_metric = squared_distance(key, item->key);
        if (item_metric <= metric) {
            memcpy(_val[num++], item->val, self->size_val);
        }

        scalar del = delta(key, item->key, depth);
        push(&top, (del < 0) ? item->left : item->right, depth + 1);
        if (sq(del) <= metric) {
            push(&top, (del < 0) ? item->right : item->left, depth + 1);
        }
    }

    arena_load(save);
    return num;
}
