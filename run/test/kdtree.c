#include "teal/kdtree.h"

#include <assert.h>
#include <stdlib.h>

#include "teal/arena.h"
#include "teal/utils.h"
#include "teal/vector.h"

static void test_basic_insert_lookup(void)
{
    Arena save = arena_save();

    Kdtree *tree = kdtree_create(sizeof(int));
    assert(tree);

    tuple num = {8, 8, 8};

    for (int i = 0; i < num.x; i++) {
        for (int j = 0; j < num.y; j++) {
            for (int k = 0; k < num.z; k++) {
                vector key = {i, j, k};
                int val = i + (num.x * (j + num.y * k));
                assert(!kdtree_insert(tree, key, &val));
            }
        }
    }

    assert(tree->num == num.x * num.y * num.z);

    for (int i = 0; i < num.x; i++) {
        for (int j = 0; j < num.y; j++) {
            for (int k = 0; k < num.z; k++) {
                vector key = {i, j, k};
                int *val = kdtree_lookup(tree, key);
                assert(val && *val == i + (num.x * (j + num.y * k)));
            }
        }
    }

    vector miss = {-1, -1, -1};
    assert(!kdtree_lookup(tree, miss));

    arena_load(save);
}

static void test_duplicate_behaviour(void)
{
    Arena save = arena_save();

    Kdtree *tree = kdtree_create(sizeof(int));
    assert(tree);

    vector key = {42, 43, 44};
    int val1 = 1111;
    int val2 = 2222;

    assert(!kdtree_insert(tree, key, &val1));
    int num = tree->num;

    int *insert = kdtree_insert(tree, key, &val2);
    assert(insert && *insert == 1111);

    assert(tree->num == num);

    int *lookup = kdtree_lookup(tree, key);
    assert(lookup == insert && *lookup == 1111);

    arena_load(save);
}

static void test_set_mode(void)
{
    Arena save = arena_save();

    Kdtree *set = kdtree_create(0);
    assert(set);

    vector key = {42, 43, 44};
    assert(!kdtree_insert(set, key, 0));
    assert(kdtree_insert(set, key, 0));
    assert(kdtree_lookup(set, key));

    vector miss = {24, 25, 26};
    assert(!kdtree_lookup(set, miss));

    arena_load(save);
}

static void test_key_by_value_not_pointer(void)
{
    Arena save = arena_save();

    Kdtree *tree = kdtree_create(sizeof(int));
    assert(tree);

    vector key = {42, 43, 44};
    int val = 1111;
    assert(!kdtree_insert(tree, key, &val));

    vector same_key = {42, 43, 44};
    int *same_val = kdtree_lookup(tree, same_key);
    assert(same_val && *same_val == 1111);

    arena_load(save);
}

static void test_pointer_stability_under_growth(void)
{
    Arena save = arena_save();

    Kdtree *tree = kdtree_create(sizeof(int));
    assert(tree);

    vector key1 = {42, 43, 44};
    int val1 = 1111;
    assert(!kdtree_insert(tree, key1, &val1));

    int *ptr1 = kdtree_lookup(tree, key1);
    assert(ptr1);

    uintptr_t addr1 = (uintptr_t)ptr1;

    for (int i = 0; i < 1000; i++) {
        vector key = {i + 420, i + 421, i + 422};
        int val = i * i;
        kdtree_insert(tree, key, &val);
    }

    int *ptr2 = kdtree_lookup(tree, key1);
    assert(ptr2);

    uintptr_t addr2 = (uintptr_t)ptr2;
    assert(addr1 == addr2);
    assert(*ptr2 == 1111);

    arena_load(save);
}

static void test_many_keys(void)
{
    Arena save = arena_save();

    Kdtree *tree = kdtree_create(sizeof(int));
    assert(tree);

    int num = 32749;                    // prime
    tuple stride = {7919, 7759, 7841};  // coprimes

    for (int i = 0; i < num; i++) {
        int x = (i * stride.x + 0) % num;
        int y = (i * stride.y + 1) % num;
        int z = (i * stride.z + 2) % num;
        vector key = {x, y, z};
        int val = x + (num * (y + num * z));
        kdtree_insert(tree, key, &val);
    }
    assert(tree->num == num);

    for (int i = 0; i < 100; i++) {
        int idx = (i * 97) % num;
        int x = (idx * stride.x + 0) % num;
        int y = (idx * stride.y + 1) % num;
        int z = (idx * stride.z + 2) % num;
        vector key = {x, y, z};
        int *val = kdtree_lookup(tree, key);
        assert(val && *val == x + (num * (y + num * z)));
    }

    arena_load(save);
}

typedef struct {
    vector key;
    int val;
    scalar dist;
} Map;

static int cmp_map(const void *lhs_, const void *rhs_)
{
    const Map *lhs = lhs_;
    const Map *rhs = rhs_;
    return cmp_asc(lhs->dist, rhs->dist);
}

static bool contains(const int *arr, int val, int num)
{
    for (int i = 0; i < num; i++) {
        if (arr[i] == val) {
            return true;
        }
    }
    return false;
}

static void test_nearest(void)
{
    Arena save = arena_save();

    Kdtree *tree = kdtree_create(sizeof(int));
    assert(tree);

    tuple cnt = {8, 8, 8};
    int tot = cnt.x * cnt.y * cnt.z;

    Map *map = arena_malloc(tot, sizeof(*map));
    int idx = 0;
    for (int i = 0; i < cnt.x; i++) {
        for (int j = 0; j < cnt.y; j++) {
            for (int k = 0; k < cnt.z; k++) {
                vector key;
                int val = i + (cnt.x * (j + cnt.y * k));
                do {
                    key = (vector){rand(), rand(), rand()};
                } while (kdtree_insert(tree, key, &val));
                map[idx].key = key;
                map[idx].val = val;
                idx += 1;
            }
        }
    }

    int num = 128;
    vector key = {rand(), rand(), rand()};
    int *val = arena_malloc(num, sizeof(*val));
    kdtree_nearest(tree, key, val, num);

    for (int i = 0; i < tot; i++) {
        vector sub = vector_sub(map[i].key, key);
        map[i].dist = vector_dot(sub, sub);
    }
    qsort(map, tot, sizeof(*map), cmp_map);

    for (int i = 0; i < num; i++) {
        assert(contains(val, map[i].val, num));
    }

    arena_load(save);
}

static void test_radius(void)
{
    Arena save = arena_save();

    Kdtree *tree = kdtree_create(sizeof(int));
    assert(tree);

    tuple cnt = {8, 8, 8};
    int tot = cnt.x * cnt.y * cnt.z;

    Map *map = arena_malloc(tot, sizeof(*map));
    int idx = 0;
    for (int i = 0; i < cnt.x; i++) {
        for (int j = 0; j < cnt.y; j++) {
            for (int k = 0; k < cnt.z; k++) {
                vector key;
                int val = i + (cnt.x * (j + cnt.y * k));
                do {
                    key = (vector){rand(), rand(), rand()};
                } while (kdtree_insert(tree, key, &val));
                map[idx].key = key;
                map[idx].val = val;
                idx += 1;
            }
        }
    }

    int cap = 128;
    scalar radius = 10;
    vector key = {rand(), rand(), rand()};
    int *val = arena_malloc(cap, sizeof(*val));
    int num = kdtree_radius(tree, key, val, cap, radius);

    for (int i = 0; i < tot; i++) {
        vector sub = vector_sub(map[i].key, key);
        map[i].dist = vector_dot(sub, sub);
    }
    qsort(map, tot, sizeof(*map), cmp_map);

    for (int i = 0; i < num; i++) {
        assert(contains(val, map[i].val, num));
    }

    arena_load(save);
}

int main(void)
{
    arena_init(str_to_size("1G"));

    test_basic_insert_lookup();
    test_duplicate_behaviour();
    test_set_mode();
    test_key_by_value_not_pointer();
    test_pointer_stability_under_growth();
    test_many_keys();
    test_nearest();
    test_radius();

    arena_finalize();
}
