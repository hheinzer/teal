#include "memory.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

static long sum = 0;

void memory_sum_setzero(void)
{
    sum = 0;
}

long memory_sum_get(void)
{
    return sum;
}

void *memory_calloc(long nmemb, long size)
{
    void *ptr = malloc(nmemb * size + sizeof(long));
    assert(ptr);
    sum += *(long *)ptr = nmemb * size;
    return memory_setzero((long *)ptr + 1, nmemb, size);
}

void *memory_realloc(void *ptr, long nmemb, long size)
{
    if (!ptr) return memory_calloc(nmemb, size);
    sum -= *((long *)ptr - 1);
    ptr = realloc((long *)ptr - 1, nmemb * size + sizeof(long));
    assert(ptr);
    sum += *(long *)ptr = nmemb * size;
    return (long *)ptr + 1;
}

void *memory_copy(void *dest, const void *src, long nmemb, long size)
{
    return memcpy(dest, src, nmemb * size);
}

void *memory_duplicate(const void *ptr, long nmemb, long size)
{
    return memory_copy(memory_calloc(nmemb, size), ptr, nmemb, size);
}

void *memory_setzero(void *ptr, long nmemb, long size)
{
    return memset(ptr, 0, nmemb * size);
}

long memory_size(const void *ptr)
{
    return *((long *)ptr - 1);
}

void memory_free(void *ptr)
{
    if (!*((void **)ptr)) return;
    sum -= *((long *)(*((void **)ptr)) - 1);
    free((long *)(*((void **)ptr)) - 1);
    *((void **)ptr) = 0;
}
