#pragma once

#include <stddef.h>

typedef struct {
    char *last;
    char *beg;
} Arena;

void arena_init(ptrdiff_t capacity);

void *arena_malloc(int num, ptrdiff_t size) __attribute((malloc));
void *arena_calloc(int num, ptrdiff_t size) __attribute((malloc));

void *arena_memdup(const void *ptr, int num, ptrdiff_t size) __attribute((malloc));

/* Resize only the most recent allocation; else returns `0`. */
void *arena_resize(const void *ptr, int num, ptrdiff_t size);

/* Move memory from a frame invalidated by `arena_load()` back into the arena. */
void *arena_smuggle(const void *ptr, int num, ptrdiff_t size) __attribute((malloc));

Arena arena_save(void);
void arena_load(Arena save);

ptrdiff_t arena_size(void);

void arena_finalize(void);
