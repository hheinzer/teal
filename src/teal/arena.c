#include "arena.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

enum { ALIGN = 64 };

#ifndef __has_feature
#define __has_feature(x) 0
#endif

#if __has_feature(address_sanitizer) || defined(__SANITIZE_ADDRESS__)

#include <sanitizer/asan_interface.h>

#define MAKE_REGION_NOACCESS(addr, size) __asan_poison_memory_region(addr, size)
#define MAKE_REGION_ADDRESSABLE(addr, size) __asan_unpoison_memory_region(addr, size)
#define MAKE_REGION_DEFINED(addr, size) MAKE_REGION_ADDRESSABLE(addr, size)

enum { REDZONE = 64 };

#elif defined(VALGRIND)

#include <valgrind/memcheck.h>

#define MAKE_REGION_NOACCESS(addr, size) VALGRIND_MAKE_MEM_NOACCESS(addr, size)
#define MAKE_REGION_ADDRESSABLE(addr, size) VALGRIND_MAKE_MEM_UNDEFINED(addr, size)
#define MAKE_REGION_DEFINED(addr, size) VALGRIND_MAKE_MEM_DEFINED(addr, size)

enum { REDZONE = 64 };

#else

#define MAKE_REGION_NOACCESS(addr, size) ((void)(addr), (void)(size))
#define MAKE_REGION_ADDRESSABLE(addr, size) ((void)(addr), (void)(size))
#define MAKE_REGION_DEFINED(addr, size) ((void)(addr), (void)(size))

enum { REDZONE = 0 };

#endif

static Arena arena = {0};
static char *base = 0;
static char *end = 0;
static long peak = 0;

void arena_init(long capacity)
{
    assert(capacity > 0);
    base = malloc(capacity);
    if (!base) {
        error("malloc failure (%ld)", capacity);
    }
    arena.beg = base;
    end = base + capacity;
    MAKE_REGION_NOACCESS(base, capacity);
}

void arena_deinit(void)
{
    free(base);
}

void *arena_malloc(long num, long size)
{
    assert(num >= 0 && size > 0);
    long available = end - arena.beg;
    long padding = -(long)arena.beg & (ALIGN - 1);
    if (num > (available - padding - REDZONE) / size) {
        error("out of memory (%ld, %ld)", num, size);
    }
    arena.last = arena.beg + padding + REDZONE;
    arena.beg = arena.last + (num * size);
    MAKE_REGION_ADDRESSABLE(arena.last, num * size);
    peak = lmax(peak, arena_size());
    return arena.last;
}

void *arena_calloc(long num, long size)
{
    assert(num >= 0 && size > 0);
    void *ptr = arena_malloc(num, size);
    return memset(ptr, 0, num * size);
}

void *arena_memdup(const void *ptr, long num, long size)
{
    assert(ptr && num >= 0 && size > 0);
    void *dup = arena_malloc(num, size);
    return memcpy(dup, ptr, num * size);
}

void *arena_resize(const void *ptr, long num, long size)
{
    if (ptr == arena.last) {
        assert(num >= 0 && size > 0);
        long available = end - arena.last;
        if (num > available / size) {
            error("out of memory (%ld, %ld)", num, size);
        }
        long old_size = arena.beg - arena.last;
        long new_size = num * size;
        arena.beg = arena.last + new_size;
        if (old_size < new_size) {
            MAKE_REGION_ADDRESSABLE(arena.last + old_size, new_size - old_size);
        }
        else {
            MAKE_REGION_NOACCESS(arena.last + new_size, old_size - new_size);
        }
        peak = lmax(peak, arena_size());
        return arena.last;
    }
    return 0;
}

static void *(*const volatile force_memmove)(void *, const void *, size_t) = memmove;

void *arena_smuggle(const void *ptr, long num, long size)
{
    assert(num >= 0 && size > 0);
    char *new = arena_malloc(num, size);
    assert(new <= (const char *)ptr && (const char *)ptr < end);
    MAKE_REGION_DEFINED(ptr, num * size);
    if (new == ptr) {
        return new;
    }
    if ((const char *)ptr - new < num * size) {
        force_memmove(new, ptr, num * size);
        MAKE_REGION_NOACCESS(new + (num * size), (const char *)ptr - new);
    }
    else {
        memcpy(new, ptr, num * size);
        MAKE_REGION_NOACCESS(ptr, num * size);
    }
    return new;
}

Arena arena_save(void)
{
    return arena;
}

void arena_load(Arena save)
{
    assert(save.beg <= arena.beg);
    MAKE_REGION_NOACCESS(save.beg, arena.beg - save.beg);
    arena = save;
}

long arena_size(void)
{
    return arena.beg - base;
}

long arena_peak(void)
{
    return peak;
}
