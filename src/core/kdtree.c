#include <assert.h>
#include <string.h>

#include "arena2.h"
#include "kdtree2.h"
#include "teal2.h"
#include "utils2.h"

enum { DIM = 3 };

typedef struct item {
    double pos[DIM];
    void *data;
    struct item *left;
    struct item *right;
} Item;

struct kdtree {
    int size;
    Item *root;
    Arena2 *arena;
};

Kdtree2 *kdtree2_init(int size)
{
    assert(size > 0);
    Kdtree2 *self = teal2_calloc(1, sizeof(*self));
    self->size = size;
    self->arena = arena2_init(0);
    return self;
}

void kdtree2_deinit(Kdtree2 *self)
{
    assert(self);
    arena2_deinit(self->arena);
    teal2_free(self);
}

static void insert_r(Kdtree2 *self, Item **item, const double *pos, const void *data, int depth)
{
    if (!*item) {
        *item = arena2_calloc(self->arena, 1, sizeof(**item));
        memcpy((*item)->pos, pos, sizeof((*item)->pos));
        (*item)->data = arena2_calloc(self->arena, 1, self->size);
        memcpy((*item)->data, data, self->size);
        return;
    }

    int dim = depth % DIM;
    item = (pos[dim] < (*item)->pos[dim]) ? &(*item)->left : &(*item)->right;
    insert_r(self, item, pos, data, depth + 1);
}

void kdtree2_insert(Kdtree2 *self, Vector pos, const void *data)
{
    assert(self && data);
    insert_r(self, &self->root, (double[]){pos.x, pos.y, pos.z}, data, 0);
}

static double distance2(const double *lhs, const double *rhs)
{
    double dist2 = 0;
    for (int i = 0; i < DIM; i++) {
        dist2 += sq(lhs[i] - rhs[i]);
    }
    return dist2;
}

typedef struct {
    const Item *item;
    double dist2;
} Hit;

static void hit_push(Hit *hit, int *num, int cap, const Item *item, double dist2)
{
    if (*num < cap) {
        int idx = *num;
        while (idx > 0 && hit[idx - 1].dist2 < dist2) {
            hit[idx] = hit[idx - 1];
            idx -= 1;
        }

        hit[idx].item = item;
        hit[idx].dist2 = dist2;
        *num += 1;
        return;
    }

    hit[0].item = item;
    hit[0].dist2 = dist2;

    int idx = 0;
    while (idx + 1 < cap && hit[idx].dist2 < hit[idx + 1].dist2) {
        Hit swap = hit[idx];
        hit[idx] = hit[idx + 1];
        hit[idx + 1] = swap;
        idx += 1;
    }
}

static void nearest_r(const Item *item, const double *pos, Hit *hit, int *num, int cap, int depth)
{
    if (!item) {
        return;
    }

    double dist2 = distance2(item->pos, pos);
    if (*num < cap || dist2 < hit[0].dist2) {
        hit_push(hit, num, cap, item, dist2);
    }

    int dim = depth % DIM;
    double delta = pos[dim] - item->pos[dim];

    Item *near = (delta < 0) ? item->left : item->right;
    nearest_r(near, pos, hit, num, cap, depth + 1);

    if (*num < cap || sq(delta) < hit[0].dist2) {
        Item *far = (delta < 0) ? item->right : item->left;
        nearest_r(far, pos, hit, num, cap, depth + 1);
    }
}

int kdtree2_nearest(const Kdtree2 *self, Vector pos, void *data_, int cap)
{
    assert(self && data_ && cap > 0);

    if (!self->root) {
        return 0;
    }

    Save2 *save = arena2_save(self->arena);
    Hit *hit = arena2_calloc(self->arena, cap, sizeof(*hit));

    int num = 0;
    nearest_r(self->root, (double[]){pos.x, pos.y, pos.z}, hit, &num, cap, 0);

    char (*data)[self->size] = data_;
    for (int i = 0; i < num; i++) {
        memcpy(data[i], hit[num - 1 - i].item->data, self->size);
    }

    arena2_load(self->arena, save);
    return num;
}

static int radius_r(const Kdtree2 *self, const Item *item, const double *pos, void *data_, int num,
                    int cap, int depth, double radius2)
{
    if (!item) {
        return num;
    }

    double dist2 = distance2(item->pos, pos);
    if (dist2 <= radius2) {
        if (num < cap) {
            char (*data)[self->size] = data_;
            memcpy(data[num], item->data, self->size);
        }
        num += 1;
    }

    int dim = depth % DIM;
    double delta = pos[dim] - item->pos[dim];

    Item *near = (delta < 0) ? item->left : item->right;
    num = radius_r(self, near, pos, data_, num, cap, depth + 1, radius2);

    if (sq(delta) <= radius2) {
        Item *far = (delta < 0) ? item->right : item->left;
        num = radius_r(self, far, pos, data_, num, cap, depth + 1, radius2);
    }

    return num;
}

int kdtree2_radius(const Kdtree2 *self, Vector pos, void *data, int cap, double radius)
{
    assert(self && data && cap > 0 && radius >= 0);
    if (!self->root) {
        return 0;
    }
    return radius_r(self, self->root, (double[]){pos.x, pos.y, pos.z}, data, 0, cap, 0, sq(radius));
}
