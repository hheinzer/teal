#pragma once

#include "teal.h"

typedef struct {
    char *last;
    char *beg;
} Arena;

void arena_init(number capacity);

void *arena_malloc(number num, number size) __attribute((malloc));
void *arena_calloc(number num, number size) __attribute((malloc));

void *arena_memdup(const void *ptr, number num, number size) __attribute((malloc));

/* Duplicate a string of length `len`; call `strlen()` if `len < 0`. */
char *arena_strdup(const char *str, number len) __attribute((malloc));

/* Resize only the most recent allocation; else returns `0`. */
void *arena_resize(const void *ptr, number num, number size);

/* Move memory from a frame invalidated by `arena_load()` back into the arena. */
void *arena_smuggle(const void *ptr, number num, number size) __attribute((malloc));

Arena arena_save(void);
void arena_load(Arena save);

number arena_size(void);
number arena_size_max(void);

void arena_finalize(void);
