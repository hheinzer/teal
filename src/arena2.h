#pragma once

#include <stddef.h>

typedef struct Arena2 {
    char *base;
    char *beg;
    char *end;
    struct Arena2 *prev;
} Arena2;

typedef struct {
    char *base;
    char *beg;
} Save2;

// Initialize an arena.
Arena2 *arena2_init(size_t capacity);

// Free an arena.
void arena2_deinit(Arena2 *self);

// Allocate `num * size` bytes from an arena and zero-initialize.
void *arena2_calloc(Arena2 *self, int num, int size) __attribute__((malloc));

// Save the current arena position for later rollback.
Save2 arena2_save(const Arena2 *self);

// Roll back the arena to a previous position.
void arena2_load(Arena2 *self, Save2 save);
