#include "memory.h"

#include <assert.h>
#include <stdlib.h>

void *memory_calloc(const long nmemb, const long size) {
    void *ptr = calloc(nmemb, size);
    assert(ptr && "calloc failure");
    return ptr;
}

void *memory_realloc(void *ptr, const long nmemb, const long size) {
    ptr = realloc(ptr, nmemb * size);
    assert(ptr && "realloc failure");
    return ptr;
}

void memory_cleanup(void *ptr) {
    free(*((void **)ptr));
}
