#pragma once

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

/* Create dictionary with enough space for 'max_items'. */
Dict dict_create(long max_items);

/* Deallocate memory of 'dict' and set all fields to 0. */
void dict_free(Dict *dict);

/* Insert 'key'-'val' pair of lengths 'nkey' and 'nval' into 'dict'. Overwrite on collision. */
void dict_insert(Dict *dict, const long *key, long nkey, const long *val, long nval);

/* Append 'key'-'val' pair of lengths 'nkey' and 'nval' in 'dict'. Insert if not present. */
void dict_append(Dict *dict, const long *key, long nkey, const long *val, long nval);

/* Lookup 'key' of length 'nkey' in 'dict', return item if present, otherwise return 0. */
DictItem *dict_lookup(const Dict *dict, const long *key, long nkey);

/* Serialize 'dict', sort items by their keys. */
DictItem *dict_serialize_by_key(const Dict *dict);

/* Serialize 'dict', sort items by their index. */
DictItem *dict_serialize_by_index(const Dict *dict);

/* Print content of 'dict'. */
void dict_print(const Dict *dict);
