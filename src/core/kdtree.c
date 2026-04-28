#include "kdtree.h"

#include <assert.h>
#include <math.h>

#include "teal.h"
#include "utils.h"

enum { CAP_LEAF = 32 };

typedef struct node {
    int beg, end;
    Vector min, max;
    struct node *left;
    struct node *right;
} Node;

struct kdtree {
    int *perm;
    Node *node;
    const Vector *point;
};

static int count_r(int num_points)
{
    if (num_points <= CAP_LEAF) {
        return 1;
    }
    int num_left = num_points / 2;
    int num_right = num_points - num_left;
    return 1 + count_r(num_left) + count_r(num_right);
}

static void compute_bbox(const Kdtree *self, Node *node)
{
    Vector min = self->point[self->perm[node->beg]];
    Vector max = min;
    for (int i = node->beg + 1; i < node->end; i++) {
        Vector point = self->point[self->perm[i]];
        min.x = fmin(min.x, point.x);
        min.y = fmin(min.y, point.y);
        min.z = fmin(min.z, point.z);
        max.x = fmax(max.x, point.x);
        max.y = fmax(max.y, point.y);
        max.z = fmax(max.z, point.z);
    }
    node->min = min;
    node->max = max;
}

static int choose_axis(const Node *node)
{
    Vector del = vector_sub(node->max, node->min);
    if (del.x >= del.y && del.x >= del.z) {
        return 0;
    }
    if (del.y >= del.z) {
        return 1;
    }
    return 2;
}

static double component(Vector vec, int axis)
{
    switch (axis) {
        case 0: return vec.x;
        case 1: return vec.y;
        case 2: return vec.z;
        default: teal_error("invalid axis (%d)", axis);
    }
}

static int compare_point(const Kdtree *self, int lhs, int rhs, int axis)
{
    Vector point_l = self->point[lhs];
    Vector point_r = self->point[rhs];
    double component_l = component(point_l, axis);
    double component_r = component(point_r, axis);
    if (component_l != component_r) {
        return (component_l < component_r) ? -1 : +1;
    }
    if (point_l.x != point_r.x) {
        return (point_l.x < point_r.x) ? -1 : +1;
    }
    if (point_l.y != point_r.y) {
        return (point_l.y < point_r.y) ? -1 : +1;
    }
    if (point_l.z != point_r.z) {
        return (point_l.z < point_r.z) ? -1 : +1;
    }
    return (lhs > rhs) - (lhs < rhs);
}

static int choose_pivot(const Kdtree *self, int beg, int end, int axis)
{
    int mid = (beg + end) / 2;
    int last = end - 1;
    int idx_beg = self->perm[beg];
    int idx_mid = self->perm[mid];
    int idx_last = self->perm[last];
    if (compare_point(self, idx_beg, idx_mid, axis) < 0) {
        if (compare_point(self, idx_mid, idx_last, axis) < 0) {
            return mid;
        }
        return (compare_point(self, idx_beg, idx_last, axis) < 0) ? last : beg;
    }
    if (compare_point(self, idx_beg, idx_last, axis) < 0) {
        return beg;
    }
    return (compare_point(self, idx_mid, idx_last, axis) < 0) ? last : mid;
}

static void swap_int(int *lhs, int *rhs)
{
    int swap = *lhs;
    *lhs = *rhs;
    *rhs = swap;
}

static int partition(Kdtree *self, int beg, int end, int axis, int pivot)
{
    int idx_pivot = self->perm[pivot];
    swap_int(&self->perm[pivot], &self->perm[end - 1]);

    int store = beg;
    for (int i = beg; i < end - 1; i++) {
        if (compare_point(self, self->perm[i], idx_pivot, axis) < 0) {
            swap_int(&self->perm[store++], &self->perm[i]);
        }
    }

    swap_int(&self->perm[store], &self->perm[end - 1]);

    return store;
}

static void quickselect(Kdtree *self, int beg, int end, int mid, int axis)
{
    while (beg + 1 < end) {
        int pivot = choose_pivot(self, beg, end, axis);
        int pos = partition(self, beg, end, axis, pivot);
        if (pos == mid) {
            return;
        }
        if (pos > mid) {
            end = pos;
        }
        else {
            beg = pos + 1;
        }
    }
}

static Node *build_nodes(Kdtree *self, Node **next, int beg, int end)
{
    Node *node = (*next)++;
    node->beg = beg;
    node->end = end;
    compute_bbox(self, node);

    if (end - beg <= CAP_LEAF) {
        return node;
    }

    int mid = (beg + end) / 2;
    int axis = choose_axis(node);
    quickselect(self, beg, end, mid, axis);

    node->left = build_nodes(self, next, beg, mid);
    node->right = build_nodes(self, next, mid, end);

    return node;
}

Kdtree *kdtree_init(const Vector *point, int num)
{
    assert((point || num == 0) && num >= 0);

    Kdtree *self = teal_calloc(1, sizeof(*self));

    if (num == 0) {
        return self;
    }

    self->point = point;
    self->perm = teal_calloc(num, sizeof(*self->perm));
    for (int i = 0; i < num; i++) {
        self->perm[i] = i;
    }

    int num_nodes = count_r(num);
    self->node = teal_calloc(num_nodes, sizeof(*self->node));

    Node *next = self->node;
    build_nodes(self, &next, 0, num);
    assert(next == self->node + num_nodes);

    return self;
}

void kdtree_deinit(Kdtree *self)
{
    assert(self);
    teal_free(self->perm);
    teal_free(self->node);
    teal_free(self);
}

static double dist2_point(Vector lhs, Vector rhs)
{
    return vector_norm2(vector_sub(lhs, rhs));
}

static double dist2_bbox(const Node *node, Vector query)
{
    double dist2 = 0;
    if (query.x < node->min.x) {
        dist2 += sq(node->min.x - query.x);
    }
    else if (node->max.x < query.x) {
        dist2 += sq(query.x - node->max.x);
    }
    if (query.y < node->min.y) {
        dist2 += sq(node->min.y - query.y);
    }
    else if (node->max.y < query.y) {
        dist2 += sq(query.y - node->max.y);
    }
    if (query.z < node->min.z) {
        dist2 += sq(node->min.z - query.z);
    }
    else if (node->max.z < query.z) {
        dist2 += sq(query.z - node->max.z);
    }
    return dist2;
}

typedef struct {
    int idx;
    double dist2;
} Hit;

static void hit_push(Hit *hit, int *num, int cap, int idx, double dist2)
{
    if (*num < cap) {
        int pos = *num;
        while (pos > 0 && hit[pos - 1].dist2 < dist2) {
            hit[pos] = hit[pos - 1];
            pos -= 1;
        }
        hit[pos].idx = idx;
        hit[pos].dist2 = dist2;
        *num += 1;
        return;
    }

    hit[0].idx = idx;
    hit[0].dist2 = dist2;

    int pos = 0;
    while (pos + 1 < cap && hit[pos].dist2 < hit[pos + 1].dist2) {
        Hit swap = hit[pos];
        hit[pos] = hit[pos + 1];
        hit[pos + 1] = swap;
        pos += 1;
    }
}

static void nearest_r(const Kdtree *self, const Node *node, Vector point, Hit *hit, int *num,
                      int cap)
{
    if (!node->left) {
        for (int i = node->beg; i < node->end; i++) {
            int idx = self->perm[i];
            double dist2 = dist2_point(self->point[idx], point);
            if (*num < cap || dist2 < hit[0].dist2) {
                hit_push(hit, num, cap, idx, dist2);
            }
        }
        return;
    }

    const Node *left = node->left;
    const Node *right = node->right;
    double dist2_left = dist2_bbox(left, point);
    double dist2_right = dist2_bbox(right, point);

    const Node *near = left;
    const Node *far = right;
    double dist2_far = dist2_right;
    if (dist2_left > dist2_right) {
        near = right;
        far = left;
        dist2_far = dist2_left;
    }

    nearest_r(self, near, point, hit, num, cap);
    if (*num < cap || dist2_far < hit[0].dist2) {
        nearest_r(self, far, point, hit, num, cap);
    }
}

int kdtree_nearest(const Kdtree *self, Vector point, int *index, int cap)
{
    assert(self && index && cap > 0);

    if (!self->node) {
        return 0;
    }

    Hit *hit = teal_calloc(cap, sizeof(*hit));

    int num = 0;
    nearest_r(self, self->node, point, hit, &num, cap);

    for (int i = 0; i < num; i++) {
        index[i] = hit[num - 1 - i].idx;
    }

    teal_free(hit);
    return num;
}
