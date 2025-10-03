#include "arena.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "assert.h"
#include "utils.h"

enum { ALIGN = 64 };

#if __has_feature(address_sanitizer) || defined(__SANITIZE_ADDRESS__)
#include <sanitizer/asan_interface.h>
#define MAKE_REGION_NOACCESS(addr, size) __asan_poison_memory_region(addr, size)
#define MAKE_REGION_ADDRESSABLE(addr, size) __asan_unpoison_memory_region(addr, size)
#define MAKE_REGION_DEFINED(addr, size) MAKE_REGION_ADDRESSABLE(addr, size)
enum { REDZONE = 64 };
#elif defined(ENABLE_VALGRIND)
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

STATIC_ASSERT(REDZONE == 0 || REDZONE % ALIGN == 0);

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
    MAKE_REGION_NOACCESS(arena_base, capacity);
}

void *arena_malloc(long num, long size)
{
    assert(num >= 0 && size > 0);

    long available = arena_end - arena.beg;
    long padding = -(uintptr_t)arena.beg & (ALIGN - 1);  // round up to next aligned location
    assert(num <= (available - padding - REDZONE) / size);

    arena.last = arena.beg + padding + REDZONE;
    arena.beg = arena.last + num * size;
    MAKE_REGION_ADDRESSABLE(arena.last, num * size);

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
            MAKE_REGION_ADDRESSABLE(arena.last + old_size, new_size - old_size);
        }
        else {
            MAKE_REGION_NOACCESS(arena.last + new_size, old_size - new_size);
        }

        size_max = lmax(size_max, arena_size());
        return arena.last;
    }
    return 0;
}

static void *(*const volatile force_memmove)(void *, const void *, size_t) = memmove;

void *arena_smuggle(const void *ptr, long num, long size)
{
    void *new = arena_malloc(num, size);
    assert((char *)new <= (char *)ptr && (char *)ptr < arena_end);
    MAKE_REGION_DEFINED(ptr, num * size);
    if (new == ptr) {
        return new;
    }
    if ((char *)ptr - (char *)new < num * size) {
        force_memmove(new, ptr, num * size);  // avoid compiler folding into memcpy
        MAKE_REGION_NOACCESS((char *)new + (num * size), (char *)ptr - (char *)new);
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
