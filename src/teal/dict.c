#include "dict.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "hash.h"
#include "memory.h"

#define LOAD_FACTOR 0.5

static int itemcmp(const void *a, const void *b);

Dict *dict_create(long max_items)
{
    Dict *dict = memory_calloc(1, sizeof(*dict));
    dict->n_items = 0;
    dict->max_items = ceil(max_items / LOAD_FACTOR) + 1;
    dict->max_dist = 1;
    dict->item = memory_calloc(dict->max_items, sizeof(*dict->item));
    return dict;
}

void dict_insert(Dict *dict, const long *key, const long *val, long nkey, long nval)
{
    assert(dict->n_items < dict->max_items * LOAD_FACTOR);
    const uint64_t hash = hash_fnv_64a(key, nkey, sizeof(*key));
    long dist = 1;
    long i = hash % dict->max_items;
    for (long n = 0; n < dict->max_items; ++n) {
        if (!dict->item[i].nkey) {
            dict->item[i].nkey = nkey;
            dict->item[i].key = memory_duplicate(key, nkey, sizeof(*key));
            dict->item[i].nval = nval;
            dict->item[i].val = memory_duplicate(val, nval, sizeof(*val));
            dict->item[i].hash = hash;
            dict->item[i].index = dict->n_items++;
            if (dist > dict->max_dist) dict->max_dist = dist;
            return;
        }
        if (dict->item[i].hash == hash) {
            memory_free(&dict->item[i].val);
            dict->item[i].nval = nval;
            dict->item[i].val = memory_duplicate(val, nval, sizeof(*val));
            return;
        }
        dist += 1;
        i = (i + 1) % dict->max_items;
    }
    abort();
}

void dict_append(Dict *dict, const long *key, const long *val, long nkey, long nval)
{
    assert(dict->n_items < dict->max_items * LOAD_FACTOR);
    const uint64_t hash = hash_fnv_64a(key, nkey, sizeof(*key));
    long dist = 1;
    long i = hash % dict->max_items;
    for (long n = 0; n < dict->max_items; ++n) {
        if (!dict->item[i].nkey) {
            dict->item[i].nkey = nkey;
            dict->item[i].key = memory_duplicate(key, nkey, sizeof(*key));
            dict->item[i].nval = nval;
            dict->item[i].val = memory_duplicate(val, nval, sizeof(*val));
            dict->item[i].hash = hash;
            dict->item[i].index = dict->n_items++;
            if (dist > dict->max_dist) dict->max_dist = dist;
            return;
        }
        if (dict->item[i].hash == hash) {
            dict->item[i].val =
                memory_realloc(dict->item[i].val, (dict->item[i].nval + nval), sizeof(*val));
            for (long j = 0; j < nval; ++j) dict->item[i].val[dict->item[i].nval + j] = val[j];
            dict->item[i].nval += nval;
            return;
        }
        dist += 1;
        i = (i + 1) % dict->max_items;
    }
    abort();
}

DictItem *dict_lookup(const Dict *dict, const long *key, long nkey)
{
    const uint64_t hash = hash_fnv_64a(key, nkey, sizeof(*key));
    long i = hash % dict->max_items;
    for (long n = 0; n < dict->max_dist; ++n) {
        if (!dict->item[i].nkey) return 0;  // WARNING: only works if items are never removed
        if (dict->item[i].hash == hash) return &dict->item[i];
        i = (i + 1) % dict->max_items;
    }
    return 0;
}

DictItem *dict_serialize_by_key(const Dict *dict)
{
    DictItem *item = memory_calloc(dict->n_items, sizeof(*item));
    for (long n = 0, i = 0; i < dict->max_items; ++i)
        if (dict->item[i].nkey) item[n++] = dict->item[i];
    qsort(item, dict->n_items, sizeof(*item), itemcmp);
    return item;
}

DictItem *dict_serialize_by_index(const Dict *dict)
{
    DictItem *item = memory_calloc(dict->n_items, sizeof(*item));
    for (long i = 0; i < dict->max_items; ++i)
        if (dict->item[i].nkey) item[dict->item[i].index] = dict->item[i];
    return item;
}

void dict_print(const Dict *dict)
{
    smart const DictItem *item = dict_serialize_by_key(dict);
    for (long i = 0; i < dict->n_items; ++i) {
        for (long j = 0; j < item[i].nkey; ++j) printf("%ld ", item[i].key[j]);
        printf(": ");
        for (long j = 0; j < item[i].nval; ++j) printf("%ld ", item[i].val[j]);
        printf("\n");
    }
}

void dict_free(Dict **dict)
{
    for (long i = 0; i < (*dict)->max_items; ++i) {
        memory_free(&(*dict)->item[i].key);
        memory_free(&(*dict)->item[i].val);
    }
    memory_free(&(*dict)->item);
    memory_free(dict);
}

static int itemcmp(const void *a, const void *b)
{
    const DictItem *aa = a;
    const DictItem *bb = b;
    if (aa->nkey < bb->nkey) return -1;
    if (aa->nkey > bb->nkey) return 1;
    for (long i = 0; i < aa->nkey; ++i) {
        if (aa->key[i] < bb->key[i]) return -1;
        if (aa->key[i] > bb->key[i]) return 1;
    }
    return 0;
}
