#include "memory.h"

#include <stdlib.h>
#include <string.h>

#include "utils.h"

void *memory_calloc(long nmemb, long size)
{
    void *ptr = calloc(nmemb, size);
    ensure(ptr);
    return ptr;
}

void *memory_realloc(void *ptr, long nmemb, long size)
{
    ptr = realloc(ptr, nmemb * size);
    ensure(ptr);
    return ptr;
}

void memory_cleanup(void *ptr)
{
    free(*((void **)ptr));
    *((void **)ptr) = 0;
}

void *memory_duplicate(const void *ptr, long nmemb, long size)
{
    return memcpy(memory_calloc(nmemb, size), ptr, nmemb * size);
}

char *memory_strdup(const char *src)
{
    return memory_duplicate(src, strlen(src) + 1, sizeof(*src));
}

void *memory_copy(void *dest, const void *src, long nmemb, long size)
{
    return memcpy(dest, src, nmemb * size);
}

void *memory_setzero(void *ptr, long nmemb, long size)
{
    return memset(ptr, 0, nmemb * size);
}
