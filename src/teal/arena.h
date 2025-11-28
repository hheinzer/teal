#pragma once

typedef struct {
    char *last;
    char *beg;
} Arena;

void arena_init(long capacity);
void arena_deinit(void);

void *arena_malloc(long num, long size) __attribute((malloc));
void *arena_calloc(long num, long size) __attribute((malloc));

void *arena_memdup(const void *ptr, long num, long size) __attribute((malloc));

// Resize only the most recent allocation; else returns `0`.
void *arena_resize(const void *ptr, long num, long size);

// Move memory from a frame invalidated by `arena_load()` back into the arena.
void *arena_smuggle(const void *ptr, long num, long size) __attribute((malloc));

Arena arena_save(void);
void arena_load(Arena save);

long arena_size(void);
long arena_peak(void);
