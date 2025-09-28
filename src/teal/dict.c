#include "dict.h"

#include <stdint.h>
#include <string.h>

#include "arena.h"
#include "utils.h"

enum { WIDTH = countof((DictItem){0}.child) };
_Static_assert(WIDTH == 2 || WIDTH == 4 || WIDTH == 8, "unsupport child array size");

enum { SHIFT = (WIDTH == 2) ? 1 : (WIDTH == 4) ? 2 : (WIDTH == 8) ? 3 : -1 };
enum { SELECT = (8 * sizeof(uint64_t)) - SHIFT };

Dict *dict_create(long size_key, long size_val)
{
    assert(size_key > 0 && size_val >= 0);
    Dict *dict = arena_calloc(1, sizeof(*dict));
    dict->end = &dict->beg;
    dict->size_key = size_key;
    dict->size_val = size_val;
    return dict;
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
