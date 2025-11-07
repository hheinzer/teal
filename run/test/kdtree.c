#include "teal/kdtree.h"

#include <stdlib.h>

#include "teal/arena.h"
#include "teal/assert.h"
#include "teal/utils.h"
#include "teal/vector.h"

static void test_basic_insert_lookup(void)
{
    Arena save = arena_save();

    Kdtree *tree = kdtree_create(sizeof(number));
    assert(tree);

    tuple num = {8, 8, 8};

    for (number i = 0; i < num.x; i++) {
        for (number j = 0; j < num.y; j++) {
            for (number k = 0; k < num.z; k++) {
                vector key = {i, j, k};
                number val = i + (num.x * (j + num.y * k));
                assert(!kdtree_insert(tree, key, &val));
            }
        }
    }

    assert(tree->num == num.x * num.y * num.z);

    for (number i = 0; i < num.x; i++) {
        for (number j = 0; j < num.y; j++) {
            for (number k = 0; k < num.z; k++) {
                vector key = {i, j, k};
                number *val = kdtree_lookup(tree, key);
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

    Kdtree *tree = kdtree_create(sizeof(number));
    assert(tree);

    vector key = {42, 43, 44};
    number val1 = 1111;
    number val2 = 2222;

    assert(!kdtree_insert(tree, key, &val1));
    number num = tree->num;

    number *insert = kdtree_insert(tree, key, &val2);
    assert(insert && *insert == 1111);

    assert(tree->num == num);

    number *lookup = kdtree_lookup(tree, key);
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

    Kdtree *tree = kdtree_create(sizeof(number));
    assert(tree);

    vector key = {42, 43, 44};
    number val = 1111;
    assert(!kdtree_insert(tree, key, &val));

    vector same_key = {42, 43, 44};
    number *same_val = kdtree_lookup(tree, same_key);
    assert(same_val && *same_val == 1111);

    arena_load(save);
}

static void test_pointer_stability_under_growth(void)
{
    Arena save = arena_save();

    Kdtree *tree = kdtree_create(sizeof(number));
    assert(tree);

    vector key1 = {42, 43, 44};
    number val1 = 1111;
    assert(!kdtree_insert(tree, key1, &val1));

    number *ptr1 = kdtree_lookup(tree, key1);
    assert(ptr1);

    number addr1 = (number)ptr1;

    for (number i = 0; i < 1000; i++) {
        vector key = {i + 420, i + 421, i + 422};
        number val = i * i;
        kdtree_insert(tree, key, &val);
    }

    number *ptr2 = kdtree_lookup(tree, key1);
    assert(ptr2);

    number addr2 = (number)ptr2;
    assert(addr1 == addr2);
    assert(*ptr2 == 1111);

    arena_load(save);
}

static void test_many_keys(void)
{
    Arena save = arena_save();

    Kdtree *tree = kdtree_create(sizeof(number));
    assert(tree);

    number num = 32749;                 // prime
    tuple stride = {7919, 7759, 7841};  // coprimes

    for (number i = 0; i < num; i++) {
        number x = (i * stride.x + 0) % num;
        number y = (i * stride.y + 1) % num;
        number z = (i * stride.z + 2) % num;
        vector key = {x, y, z};
        number val = x + (num * (y + num * z));
        kdtree_insert(tree, key, &val);
    }
    assert(tree->num == num);

    for (number i = 0; i < 100; i++) {
        number idx = (i * 97) % num;
        number x = (idx * stride.x + 0) % num;
        number y = (idx * stride.y + 1) % num;
        number z = (idx * stride.z + 2) % num;
        vector key = {x, y, z};
        number *val = kdtree_lookup(tree, key);
        assert(val && *val == x + (num * (y + num * z)));
    }

    arena_load(save);
}

typedef struct {
    vector key;
    number val;
    scalar dist;
} Map;

static int cmp_map(const void *lhs_, const void *rhs_)
{
    const Map *lhs = lhs_;
    const Map *rhs = rhs_;
    return cmp_asc(lhs->dist, rhs->dist);
}

static bool contains(const number *arr, number val, number num)
{
    for (number i = 0; i < num; i++) {
        if (arr[i] == val) {
            return true;
        }
    }
    return false;
}

static void test_nearest(void)
{
    Arena save = arena_save();

    Kdtree *tree = kdtree_create(sizeof(number));
    assert(tree);

    tuple cnt = {8, 8, 8};
    number tot = cnt.x * cnt.y * cnt.z;

    Map *map = arena_malloc(tot, sizeof(*map));
    number idx = 0;
    for (number i = 0; i < cnt.x; i++) {
        for (number j = 0; j < cnt.y; j++) {
            for (number k = 0; k < cnt.z; k++) {
                vector key;
                number val = i + (cnt.x * (j + cnt.y * k));
                do {
                    key = (vector){rand(), rand(), rand()};
                } while (kdtree_insert(tree, key, &val));
                map[idx].key = key;
                map[idx].val = val;
                idx += 1;
            }
        }
    }

    number num = 128;
    vector key = {rand(), rand(), rand()};
    number *val = arena_malloc(num, sizeof(*val));
    kdtree_nearest(tree, key, val, num);

    for (number i = 0; i < tot; i++) {
        vector sub = vector_sub(map[i].key, key);
        map[i].dist = vector_dot(sub, sub);
    }
    qsort(map, tot, sizeof(*map), cmp_map);

    for (number i = 0; i < num; i++) {
        assert(contains(val, map[i].val, num));
    }

    arena_load(save);
}

static void test_radius(void)
{
    Arena save = arena_save();

    Kdtree *tree = kdtree_create(sizeof(number));
    assert(tree);

    tuple cnt = {8, 8, 8};
    number tot = cnt.x * cnt.y * cnt.z;

    Map *map = arena_malloc(tot, sizeof(*map));
    number idx = 0;
    for (number i = 0; i < cnt.x; i++) {
        for (number j = 0; j < cnt.y; j++) {
            for (number k = 0; k < cnt.z; k++) {
                vector key;
                number val = i + (cnt.x * (j + cnt.y * k));
                do {
                    key = (vector){rand(), rand(), rand()};
                } while (kdtree_insert(tree, key, &val));
                map[idx].key = key;
                map[idx].val = val;
                idx += 1;
            }
        }
    }

    number cap = 128;
    scalar radius = 10;
    vector key = {rand(), rand(), rand()};
    number *val = arena_malloc(cap, sizeof(*val));
    number num = kdtree_radius(tree, key, val, cap, radius);

    for (number i = 0; i < tot; i++) {
        vector sub = vector_sub(map[i].key, key);
        map[i].dist = vector_dot(sub, sub);
    }
    qsort(map, tot, sizeof(*map), cmp_map);

    for (number i = 0; i < num; i++) {
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
