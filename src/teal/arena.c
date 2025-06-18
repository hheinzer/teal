#include "arena.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

enum { ALIGN = 64 };  // getconf LEVEL1_DCACHE_LINESIZE

#if __has_feature(address_sanitizer) || defined(__SANITIZE_ADDRESS__)
#include <sanitizer/asan_interface.h>
#define ASAN_POISON_MEMORY_REGION(addr, size) __asan_poison_memory_region((addr), (size))
#define ASAN_UNPOISON_MEMORY_REGION(addr, size) __asan_unpoison_memory_region((addr), (size))
enum { REDZONE = ALIGN };
#else
#define ASAN_POISON_MEMORY_REGION(addr, size) ((void)(addr), (void)(size))
#define ASAN_UNPOISON_MEMORY_REGION(addr, size) ((void)(addr), (void)(size))
enum { REDZONE = 0 };
#endif

static Arena arena = {0};
static char *arena_base = 0;
static char *arena_end = 0;
static long size_max = 0;

void arena_init(long capacity)
{
    assert(capacity > 0);

    arena_base = malloc(capacity);
    assert(arena_base);

    arena.beg = arena_base;
    arena_end = arena_base + capacity;

    ASAN_POISON_MEMORY_REGION(arena_base, capacity);
}

void *arena_malloc(long num, long size)
{
    assert(num >= 0 && size > 0);

    long available = arena_end - arena.beg;
    long padding = -(uintptr_t)arena.beg & (ALIGN - 1);  // round up to next aligned location
    assert(num <= (available - padding - REDZONE) / size);

    arena.last = arena.beg + padding + REDZONE;
    arena.beg = arena.last + num * size;

    ASAN_UNPOISON_MEMORY_REGION(arena.last, num * size);

    size_max = lmax(size_max, arena_size());
    return arena.last;
}

void *arena_calloc(long num, long size)
{
    void *ptr = arena_malloc(num, size);
    return memset(ptr, 0, num * size);
}

void *arena_memdup(const void *ptr, long num, long size)
{
    assert(ptr);
    void *new = arena_malloc(num, size);
    return memcpy(new, ptr, num * size);
}

void *arena_resize(const void *ptr, long num, long size)
{
    if (ptr == arena.last) {
        assert(num >= 0 && size > 0);

        long available = arena_end - arena.last;
        assert(num <= available / size);

        long old_size = arena.beg - arena.last;
        long new_size = num * size;
        arena.beg = arena.last + new_size;

        if (old_size < new_size) {
            ASAN_UNPOISON_MEMORY_REGION(arena.last + old_size, new_size - old_size);
        }
        else {
            ASAN_POISON_MEMORY_REGION(arena.last + new_size, old_size - new_size);
        }

        size_max = lmax(size_max, arena_size());
        return arena.last;
    }
    return 0;
}

void *arena_smuggle(const void *ptr, long num, long size)
{
    assert(arena.beg <= (char *)ptr && (char *)ptr < arena_end);
    void *new = arena_malloc(num, size);
    ASAN_UNPOISON_MEMORY_REGION(ptr, num * size);
    memmove(new, ptr, num * size);
    ASAN_POISON_MEMORY_REGION(ptr, num * size);
    return new;
}

void *arena_consume(void *ptr, long num, long size)
{
    assert((char *)ptr < arena_base || arena_end <= (char *)ptr);
    void *new = arena_memdup(ptr, num, size);
    free(ptr);
    return new;
}

Arena arena_save(void)
{
    return arena;
}

void arena_load(Arena save)
{
    assert(save.beg <= arena.beg);
    ASAN_POISON_MEMORY_REGION(save.beg, arena.beg - save.beg);
    arena = save;
}

long arena_size(void)
{
    return arena.beg - arena_base;
}

long arena_size_max(void)
{
    return size_max;
}

void arena_finalize(void)
{
    free(arena_base);
}
