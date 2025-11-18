#pragma once

#include "teal.h"

typedef struct Dict Dict;
typedef struct DictItem DictItem;

struct Dict {
    int num;
    DictItem *beg;
    DictItem **end;
    int size_key;
    int size_val;
};

struct DictItem {
    void *key;
    void *val;
    DictItem *child[4];
    DictItem *next;
};

/* Returns empty dictionary; if `size_val == 0` it behaves as a set (values ignored). */
Dict *dict_create(int size_key, int size_val);

/* Returns pointer to stored value if `key` is already present; else inserts and returns `0`. */
void *dict_insert(Dict *self, const void *key, const void *val);

/* Returns pointer to stored value if `key` is present or `0` if absent. */
void *dict_lookup(const Dict *self, const void *key);
