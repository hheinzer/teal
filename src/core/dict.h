#ifndef DICT_H
#define DICT_H

#include <stdint.h>

typedef struct DictItem {
    long nkey, *key;
    long nval, *val;
    uint64_t hash;
    long index;
} DictItem;

typedef struct Dict {
    long n_items, max_items, max_dist;
    DictItem *item;
} Dict;

Dict dict_create(const long n_items);

void dict_free(Dict *dict);

void dict_insert(Dict *dict, const long *key, const long nkey, const long *val, const long nval);

void dict_append(Dict *dict, const long *key, const long nkey, const long *val, const long nval);

long dict_lookup(const Dict *dict, const long *key, const long nkey, long **val);

DictItem *dict_serialize_by_key(const Dict *dict);

DictItem *dict_serialize_by_index(const Dict *dict);

void dict_print(const Dict *dict);

#endif
