#include "teal/dict.h"

#include "teal/arena.h"
#include "teal/assert.h"
#include "teal/utils.h"

static void test_basic_insert_lookup(void)
{
    Arena save = arena_save();

    Dict *dict = dict_create(sizeof(number), sizeof(number));
    assert(dict);

    number num = 1024;

    for (number i = 0; i < num; i++) {
        number key = i;
        number val = i * 10;
        assert(!dict_insert(dict, &key, &val));
    }

    assert(dict->num == num);

    for (number i = 0; i < num; i++) {
        number key = i;
        number *val = dict_lookup(dict, &key);
        assert(val && *val == i * 10);
    }

    number miss = -1;
    assert(!dict_lookup(dict, &miss));

    arena_load(save);
}

static void test_duplicate_behaviour(void)
{
    Arena save = arena_save();

    Dict *dict = dict_create(sizeof(number), sizeof(number));
    assert(dict);

    number key = 42;
    number val1 = 1111;
    number val2 = 2222;

    assert(!dict_insert(dict, &key, &val1));
    number num = dict->num;

    number *insert = dict_insert(dict, &key, &val2);
    assert(insert && *insert == 1111);

    assert(dict->num == num);

    number *lookup = dict_lookup(dict, &key);
    assert(lookup == insert && *lookup == 1111);

    arena_load(save);
}

static void test_set_mode(void)
{
    Arena save = arena_save();

    Dict *set = dict_create(sizeof(number), 0);
    assert(set);

    number key = 42;
    assert(!dict_insert(set, &key, 0));
    assert(dict_insert(set, &key, 0));
    assert(dict_lookup(set, &key));

    number miss = 24;
    assert(!dict_lookup(set, &miss));

    arena_load(save);
}

static void test_key_by_value_not_pointer(void)
{
    Arena save = arena_save();

    Dict *dict = dict_create(sizeof(number), sizeof(number));
    assert(dict);

    number key = 42;
    number val = 1111;
    assert(!dict_insert(dict, &key, &val));

    number same_key = 42;
    number *same_val = dict_lookup(dict, &same_key);
    assert(same_val && *same_val == 1111);

    arena_load(save);
}

static void test_pointer_stability_under_growth(void)
{
    Arena save = arena_save();

    Dict *dict = dict_create(sizeof(number), sizeof(number));
    assert(dict);

    number key1 = 42;
    number val1 = 1111;
    assert(!dict_insert(dict, &key1, &val1));

    number *ptr1 = dict_lookup(dict, &key1);
    assert(ptr1);

    number addr1 = (number)ptr1;

    for (number i = 0; i < 1000; i++) {
        number key = i + 420;
        number val = i * i;
        dict_insert(dict, &key, &val);
    }

    number *ptr2 = dict_lookup(dict, &key1);
    assert(ptr2);

    number addr2 = (number)ptr2;
    assert(addr1 == addr2);
    assert(*ptr2 == 1111);

    arena_load(save);
}

static void test_many_keys(void)
{
    Arena save = arena_save();

    Dict *dict = dict_create(sizeof(number), sizeof(number));
    assert(dict);

    number num = 32749;    // prime
    number stride = 7919;  // coprime

    for (number i = 0; i < num; i++) {
        number key = (i * stride) % num;
        number val = key * 3;
        dict_insert(dict, &key, &val);
    }
    assert(dict->num == num);

    for (number i = 0; i < 100; i++) {
        number key = (i * 97) % num;
        number *val = dict_lookup(dict, &key);
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
