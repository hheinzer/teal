#pragma once

#include <stddef.h>

typedef struct arena Arena;

typedef struct save Save;

// Initialize an arena.
Arena *arena2_init(ptrdiff_t capacity);

// Free an arena.
void arena2_deinit(Arena *self);

// Allocate `num * size` bytes from an arena and zero-initialize.
void *arena2_calloc(Arena *self, int num, int size) __attribute__((malloc));

// Save the current arena position for later rollback.
Save *arena2_save(Arena *self);

// Roll back the arena to a previous position.
void arena2_load(Arena *self, const Save *save);
