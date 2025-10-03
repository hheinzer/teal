#include <stdint.h>
#include <stdlib.h>

#include "teal.h"
#include "teal/arena.h"
#include "teal/assert.h"
#include "teal/dict.h"
#include "teal/sync.h"

static void test_basic_insert_lookup(void)
{
    Arena save = arena_save();

    Dict *dict = dict_create(sizeof(long), sizeof(long));
    assert(dict);

    long num = 1024;

    for (long i = 0; i < num; i++) {
        long key = i;
        long val = i * 10;
        assert(!dict_insert(dict, &key, &val));
    }

    assert(dict->num == num);

    for (long i = 0; i < num; i++) {
        long key = i;
        long *val = dict_lookup(dict, &key);
        assert(val && *val == i * 10);
    }

    long miss = -1;
    assert(!dict_lookup(dict, &miss));

    arena_load(save);
}

static void test_duplicate_behaviour(void)
{
    Arena save = arena_save();

    Dict *dict = dict_create(sizeof(long), sizeof(long));
    assert(dict);

    long key = 42;
    long val1 = 1111;
    long val2 = 2222;

    assert(!dict_insert(dict, &key, &val1));
    long num = dict->num;

    long *insert = dict_insert(dict, &key, &val2);
    assert(insert && *insert == 1111);

    assert(dict->num == num);

    long *lookup = dict_lookup(dict, &key);
    assert(lookup == insert && *lookup == 1111);

    arena_load(save);
}

static void test_set_mode(void)
{
    Arena save = arena_save();

    Dict *set = dict_create(sizeof(long), 0);
    assert(set);

    long key = 42;
    assert(!dict_insert(set, &key, 0));
    assert(dict_insert(set, &key, 0));
    assert(dict_lookup(set, &key));

    long miss = 24;
    assert(!dict_lookup(set, &miss));

    arena_load(save);
}

static void test_key_by_value_not_pointer(void)
{
    Arena save = arena_save();

    Dict *dict = dict_create(sizeof(long), sizeof(long));
    assert(dict);

    long key = 42;
    long val = 1111;
    assert(!dict_insert(dict, &key, &val));

    long same_key = 42;
    long *same_val = dict_lookup(dict, &same_key);
    assert(same_val && *same_val == 1111);

    arena_load(save);
}

static void test_pointer_stability_under_growth(void)
{
    Arena save = arena_save();

    Dict *dict = dict_create(sizeof(long), sizeof(long));
    assert(dict);

    long key1 = 42;
    long val1 = 1111;
    assert(!dict_insert(dict, &key1, &val1));

    long *ptr1 = dict_lookup(dict, &key1);
    assert(ptr1);

    uintptr_t addr1 = (uintptr_t)ptr1;

    for (long i = 0; i < 50000; i++) {
        long key = i + 420;
        long val = i * i;
        dict_insert(dict, &key, &val);
    }

    long *ptr2 = dict_lookup(dict, &key1);
    assert(ptr2);

    uintptr_t addr2 = (uintptr_t)ptr2;
    assert(addr1 == addr2);
    assert(*ptr2 == 1111);

    arena_load(save);
}

static void test_many_keys(void)
{
    Arena save = arena_save();

    Dict *dict = dict_create(sizeof(long), sizeof(long));
    assert(dict);

    long num = 32749;    // prime
    long stride = 7919;  // coprime

    for (long i = 0; i < num; i++) {
        long key = (i * stride) % num;
        long val = key * 3;
        dict_insert(dict, &key, &val);
    }
    assert(dict->num == num);

    for (long i = 0; i < 100; i++) {
        long key = (i * 97) % num;
        long *val = dict_lookup(dict, &key);
        assert(val && *val == key * 3);
    }

    arena_load(save);
}

static long hash(long key)
{
    long val = key;
    val ^= val >> 23;
    val *= 0x9e3779b1;
    val ^= val >> 16;
    return val;
}

static void test_stress(long num)
{
    Arena save = arena_save();

    Dict *dict = dict_create(sizeof(long), sizeof(long));
    for (long i = 0; i < num; i++) {
        long key = rand();
        long val = hash(key);
        dict_insert(dict, &key, &val);
    }

    for (long i = 0; i < num; i++) {
        long key = rand();
        long *val = dict_lookup(dict, &key);
        if (val) {
            assert(*val == hash(key));
        }
    }

    arena_load(save);
}

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);

    srand(sync.rank);

    test_basic_insert_lookup();
    test_duplicate_behaviour();
    test_set_mode();
    test_key_by_value_not_pointer();
    test_pointer_stability_under_growth();
    test_many_keys();

    for (long i = 0; i < 10; i++) {
        test_stress(1 << 21);
    }

    teal_finalize();
}
