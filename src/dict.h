#ifndef DICT_H
#define DICT_H

typedef struct DictItem {
    long n_key, n_val;
    long *key, *val;
} DictItem;

typedef struct Dict {
    long n_items, max_items, max_dist;
    DictItem *item;
} Dict;

Dict dict_create(const long n_items);

void dict_free(Dict *dict);

void dict_insert(Dict *dict, const long *key, const long n_key, const long *val, const long n_val);

void dict_append(Dict *dict, const long *key, const long n_key, const long *val, const long n_val);

long dict_lookup(const Dict *dict, const long *key, const long n_key, long **val);

DictItem *dict_serialize(const Dict *dict);

void dict_print(const Dict *dict);

#endif
