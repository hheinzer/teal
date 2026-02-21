#pragma once

#include <stddef.h>

// Arena allocator with checkpoint rollback support.
typedef struct arena Arena2;

// Checkpoint token.
typedef struct save Save2;

// Initialize an arena.
Arena2 *arena2_init(size_t capacity);

// Free an arena.
void arena2_deinit(Arena2 *self);

// Allocate `num * size` bytes from an arena and zero-initialize.
void *arena2_calloc(Arena2 *self, int num, int size) __attribute__((malloc));

// Save the current arena position for later rollback.
Save2 *arena2_save(Arena2 *self);

// Roll back the arena to a previous position.
void arena2_load(Arena2 *self, const Save2 *save);
