#include "dict.h"

#include <assert.h>
#include <stdint.h>
#include <string.h>

#include "arena.h"

#define NON_NULL ((void *)(int)1)  // NOLINT(performance-no-int-to-ptr)

#define countof(arr) (sizeof(arr) / sizeof(*(arr)))

enum { WIDTH = countof((DictItem){0}.child) };
enum { SHIFT = (WIDTH == 2) ? 1 : (WIDTH == 4) ? 2 : (WIDTH == 8) ? 3 : -1 };
enum { SELECT = (8 * sizeof(uint64_t)) - SHIFT };

Dict *dict_create(int size_key, int size_val)
{
    assert(size_key > 0 && size_val >= 0);
    Dict *dict = arena_calloc(1, sizeof(*dict));
    dict->end = &dict->beg;
    dict->size_key = size_key;
    dict->size_val = size_val;
    return dict;
}

static uint64_t fnv1a(const void *ptr, int size)
{
    static const uint64_t basis = 0xcbf29ce484222325;
    static const uint64_t prime = 0x00000100000001b3;
    const char *byte = ptr;
    uint64_t hash = basis;
    for (int i = 0; i < size; i++) {
        hash ^= byte[i];
        hash *= prime;
    }
    return hash;
}

void *dict_insert(Dict *self, const void *key, const void *val)
{
    assert(self && key && (val || self->size_val == 0));

    DictItem **item = &self->beg;
    for (uint64_t hash = fnv1a(key, self->size_key); *item; hash <<= SHIFT) {
        if (!memcmp((*item)->key, key, self->size_key)) {
            return (*item)->val;
        }
        item = &(*item)->child[hash >> SELECT];
    }

    *item = arena_calloc(1, sizeof(**item));
    (*item)->key = arena_memdup(key, 1, self->size_key);
    (*item)->val = self->size_val ? arena_memdup(val, 1, self->size_val) : NON_NULL;

    self->num += 1;
    *self->end = *item;
    self->end = &(*item)->next;

    return 0;
}

void *dict_lookup(const Dict *self, const void *key)
{
    assert(self && key);
    DictItem *item = self->beg;
    for (uint64_t hash = fnv1a(key, self->size_key); item; hash <<= SHIFT) {
        if (!memcmp(item->key, key, self->size_key)) {
            return item->val;
        }
        item = item->child[hash >> SELECT];
    }
    return 0;
}
