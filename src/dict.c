#include "dict.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "array.h"
#include "hash.h"
#include "memory.h"
#include "utils.h"

#define LOAD_FACTOR 0.5

static int keycmp(const long *key1, const long n_key1, const long *key2, const long n_key2);
static int itemcmp(const void *a, const void *b);

Dict dict_create(const long n_items) {
    Dict dict = {
        .n_items = 0,
        .max_items = ceil(n_items / LOAD_FACTOR),
        .max_dist = 1,
    };
    dict.item = memory_calloc(dict.max_items, sizeof(*dict.item));
    return dict;
}

void dict_free(Dict *dict) {
    for (long i = 0; i < dict->max_items; ++i) {
        free(dict->item[i].key);
        free(dict->item[i].val);
    }
    free(dict->item);
    *dict = (typeof(*dict)){};
}

void dict_insert(Dict *dict, const long *key, const long n_key, const long *val, const long n_val) {
    assert(dict->n_items < dict->max_items && "out of capacity");
    long dist = 1;
    long i = hash_fnv_64a(key, n_key, sizeof(*key)) % dict->max_items;
    for (long n = 0; n < dict->max_items; ++n) {
        if (!dict->item[i].key) {
            // insert new value
            dict->item[i].n_key = n_key;
            dict->item[i].n_val = n_val;
            dict->item[i].key = utils_memdup(key, n_key, sizeof(*key));
            dict->item[i].val = utils_memdup(val, n_val, sizeof(*val));
            dict->n_items += 1;
            dict->max_dist = MAX(dict->max_dist, dist);
            return;
        } else if (!keycmp(dict->item[i].key, dict->item[i].n_key, key, n_key)) {
            // overwrite existing value
            free(dict->item[i].val);
            dict->item[i].n_val = n_val;
            dict->item[i].val = utils_memdup(val, n_val, sizeof(*val));
            return;
        }
        dist += 1;
        i = (i + 1) % dict->max_items;
    }
}

void dict_append(Dict *dict, const long *key, const long n_key, const long *val, const long n_val) {
    assert(dict->n_items < dict->max_items && "out of capacity");
    long dist = 1;
    long i = hash_fnv_64a(key, n_key, sizeof(*key)) % dict->max_items;
    for (long n = 0; n < dict->max_items; ++n) {
        if (!dict->item[i].key) {
            // insert new value
            dict->item[i].n_key = n_key;
            dict->item[i].n_val = n_val;
            dict->item[i].key = utils_memdup(key, n_key, sizeof(*key));
            dict->item[i].val = utils_memdup(val, n_val, sizeof(*val));
            dict->n_items += 1;
            dict->max_dist = MAX(dict->max_dist, dist);
            return;
        } else if (!keycmp(dict->item[i].key, dict->item[i].n_key, key, n_key)) {
            // append to existing value
            dict->item[i].val =
                memory_realloc(dict->item[i].val, dict->item[i].n_val + n_val, sizeof(*val));
            for (long j = 0; j < n_val; ++j) dict->item[i].val[dict->item[i].n_val + j] = val[j];
            dict->item[i].n_val += n_val;
            return;
        }
        dist += 1;
        i = (i + 1) % dict->max_items;
    }
}

long dict_lookup(const Dict *dict, const long *key, const long n_key, long **val) {
    *val = 0;
    long i = hash_fnv_64a(key, n_key, sizeof(*key)) % dict->max_items;
    for (long n = 0; n < dict->max_dist; ++n) {
        if (!dict->item[i].key) {
            return 0;  // key does not exists
        } else if (!keycmp(dict->item[i].key, dict->item[i].n_key, key, n_key)) {
            // key found
            *val = dict->item[i].val;
            return dict->item[i].n_val;
        }
        i = (i + 1) % dict->max_items;
    }
    return 0;  // key does not exists
}

DictItem *dict_serialize(const Dict *dict) {
    DictItem *item = memory_calloc(dict->n_items, sizeof(*item));
    for (long n = 0, i = 0; i < dict->max_items; ++i) {
        if (dict->item[i].key) {
            item[n].n_key = dict->item[i].n_key;
            item[n].n_val = dict->item[i].n_val;
            item[n].key = dict->item[i].key;
            item[n].val = dict->item[i].val;
            n += 1;
        }
    }
    qsort(item, dict->n_items, sizeof(*item), itemcmp);
    return item;
}

void dict_print(const Dict *dict) {
    cleanup const DictItem *item = dict_serialize(dict);
    for (long i = 0; i < dict->n_items; ++i) {
        array_print(0, item[i].key, item[i].n_key, ": ");
        array_print(0, item[i].val, item[i].n_val, "\n");
    }
}

static int keycmp(const long *key1, const long n_key1, const long *key2, const long n_key2) {
    if (n_key1 < n_key2) return -1;
    if (n_key1 > n_key2) return 1;
    for (long i = 0; i < n_key1; ++i) {
        if (key1[i] < key2[i]) return -1;
        if (key1[i] > key2[i]) return 1;
    }
    return 0;
}

static int itemcmp(const void *a, const void *b) {
    const DictItem *item_a = (typeof(item_a))a;
    const DictItem *item_b = (typeof(item_b))b;
    return keycmp(item_a->key, item_a->n_key, item_b->key, item_b->n_key);
}
