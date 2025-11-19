#include "teal/dict.h"

#include <assert.h>

#include "teal/arena.h"
#include "teal/utils.h"

static void test_basic_insert_lookup(void)
{
    Arena save = arena_save();

    Dict *dict = dict_create(sizeof(int), sizeof(int));
    assert(dict);

    int num = 1024;

    for (int i = 0; i < num; i++) {
        int key = i;
        int val = i * 10;
        assert(!dict_insert(dict, &key, &val));
    }

    assert(dict->num == num);

    for (int i = 0; i < num; i++) {
        int key = i;
        int *val = dict_lookup(dict, &key);
        assert(val && *val == i * 10);
    }

    int miss = -1;
    assert(!dict_lookup(dict, &miss));

    arena_load(save);
}

static void test_duplicate_behaviour(void)
{
    Arena save = arena_save();

    Dict *dict = dict_create(sizeof(int), sizeof(int));
    assert(dict);

    int key = 42;
    int val1 = 1111;
    int val2 = 2222;

    assert(!dict_insert(dict, &key, &val1));
    int num = dict->num;

    int *insert = dict_insert(dict, &key, &val2);
    assert(insert && *insert == 1111);

    assert(dict->num == num);

    int *lookup = dict_lookup(dict, &key);
    assert(lookup == insert && *lookup == 1111);

    arena_load(save);
}

static void test_set_mode(void)
{
    Arena save = arena_save();

    Dict *set = dict_create(sizeof(int), 0);
    assert(set);

    int key = 42;
    assert(!dict_insert(set, &key, 0));
    assert(dict_insert(set, &key, 0));
    assert(dict_lookup(set, &key));

    int miss = 24;
    assert(!dict_lookup(set, &miss));

    arena_load(save);
}

static void test_key_by_value_not_pointer(void)
{
    Arena save = arena_save();

    Dict *dict = dict_create(sizeof(int), sizeof(int));
    assert(dict);

    int key = 42;
    int val = 1111;
    assert(!dict_insert(dict, &key, &val));

    int same_key = 42;
    int *same_val = dict_lookup(dict, &same_key);
    assert(same_val && *same_val == 1111);

    arena_load(save);
}

static void test_pointer_stability_under_growth(void)
{
    Arena save = arena_save();

    Dict *dict = dict_create(sizeof(int), sizeof(int));
    assert(dict);

    int key1 = 42;
    int val1 = 1111;
    assert(!dict_insert(dict, &key1, &val1));

    int *ptr1 = dict_lookup(dict, &key1);
    assert(ptr1);

    uintptr_t addr1 = (uintptr_t)ptr1;

    for (int i = 0; i < 1000; i++) {
        int key = i + 420;
        int val = i * i;
        dict_insert(dict, &key, &val);
    }

    int *ptr2 = dict_lookup(dict, &key1);
    assert(ptr2);

    uintptr_t addr2 = (uintptr_t)ptr2;
    assert(addr1 == addr2);
    assert(*ptr2 == 1111);

    arena_load(save);
}

static void test_many_keys(void)
{
    Arena save = arena_save();

    Dict *dict = dict_create(sizeof(int), sizeof(int));
    assert(dict);

    int num = 32749;    // prime
    int stride = 7919;  // coprime

    for (int i = 0; i < num; i++) {
        int key = (i * stride) % num;
        int val = key * 3;
        dict_insert(dict, &key, &val);
    }
    assert(dict->num == num);

    for (int i = 0; i < 100; i++) {
        int key = (i * 97) % num;
        int *val = dict_lookup(dict, &key);
        assert(val && *val == key * 3);
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

    arena_finalize();
}
