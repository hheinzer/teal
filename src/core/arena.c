#include <assert.h>
#include <string.h>

#include "arena2.h"
#include "sanitizer.h"
#include "teal2.h"

enum { ALIGN = 64 };

struct Arena2 {
    char *base;
    char *beg;
    char *end;
    Arena2 *prev;
};

struct Save2 {
    char *base;
    char *beg;
};

Arena2 *arena2_init(size_t capacity)
{
    size_t min_capacity = 10 << 20;
    if (capacity < min_capacity) {
        capacity = min_capacity;
    }

    Arena2 *self = teal2_malloc(sizeof(*self) + capacity);

    self->base = (char *)(self + 1);
    MAKE_REGION_NOACCESS(self->base, capacity);

    self->beg = self->base;
    self->end = self->base + capacity;
    self->prev = 0;

    return self;
}

void arena2_deinit(Arena2 *self)
{
    while (self) {
        Arena2 *prev = self->prev;
        teal2_free(self);
        self = prev;
    }
}

static void append_chunk(Arena2 *self, int num, int size, int align)
{
    if (num > (PTRDIFF_MAX - REDZONE - (align - 1)) / size) {
        teal2_error("overflow (%d, %d)", num, size);
    }

    ptrdiff_t old_capacity = self->end - self->base;
    ptrdiff_t min_capacity = REDZONE + (align - 1) + ((ptrdiff_t)num * size);

    if (old_capacity > min_capacity) {
        min_capacity = old_capacity;
    }

    if (min_capacity > PTRDIFF_MAX / 2) {
        teal2_error("overflow (%td)", min_capacity);
    }

    ptrdiff_t capacity = 2 * min_capacity;
    teal2_verbose("growing scratch by %td MiB", capacity >> 20);

    Arena2 *next = arena2_init(capacity);

    Arena2 swap = *self;
    *self = *next;
    *next = swap;

    self->prev = next;
}

static void *arena2_malloc(Arena2 *self, int num, int size)
{
    assert(self && num >= 0 && size > 0);

    if (num == 0) {
        return 0;
    }

    ptrdiff_t available = self->end - self->beg;
    if (num > (available - REDZONE - (ALIGN - 1)) / size) {
        append_chunk(self, num, size, ALIGN);
    }

    ptrdiff_t padding = -(intptr_t)(self->beg + REDZONE) & (ALIGN - 1);
    char *ptr = self->beg + REDZONE + padding;

    ptrdiff_t bytes = (ptrdiff_t)num * size;
    self->beg = ptr + bytes;

    MAKE_REGION_ADDRESSABLE(ptr, bytes);

    return ptr;
}

void *arena2_calloc(Arena2 *self, int num, int size)
{
    assert(self && num >= 0 && size > 0);

    if (num == 0) {
        return 0;
    }

    void *ptr = arena2_malloc(self, num, size);

    return memset(ptr, 0, (size_t)num * size);
}

void *arena2_memdup(Arena2 *self, const void *ptr, int num, int size)
{
    assert(self && (ptr || num == 0) && num >= 0 && size > 0);

    if (num == 0) {
        return 0;
    }

    void *dup = arena2_malloc(self, num, size);

    return memcpy(dup, ptr, (size_t)num * size);
}

Save2 *arena2_save(Arena2 *self)
{
    assert(self);

    char *base = self->base;
    char *beg = self->beg;

    Save2 *save = arena2_calloc(self, 1, sizeof(*save));

    save->base = base;
    save->beg = beg;

    return save;
}

void arena2_load(Arena2 *self, const Save2 *save)
{
    assert(self && save);

    while (self->prev && self->base != save->base) {
        Arena2 *prev = self->prev;
        *self = *prev;
        teal2_free(prev);
    }

    int matching_base = (self->base == save->base);
    int inside_base = (self->base <= save->beg && save->beg <= self->end);
    int from_past = (save->beg <= self->beg);
    if (!matching_base || !inside_base || !from_past) {
        teal2_error("invalid checkpoint (%d, %d, %d)", matching_base, inside_base, from_past);
    }

    char *beg = save->beg;

    assert(self->beg > save->beg);
    MAKE_REGION_NOACCESS(save->beg, self->beg - save->beg);

    self->beg = beg;
}
