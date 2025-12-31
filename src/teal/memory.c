#include <assert.h>
#include <stdint.h>
#include <stdlib.h>

#include "teal.h"

void *teal_alloc(int num, size_t size)
{
    assert(num >= 0 && size > 0);
    if (num == 0) {
        return 0;
    }
    if ((size_t)num > SIZE_MAX / size) {
        teal_error("size overflow");
    }
    void *ptr = calloc((size_t)num, size);
    if (!ptr) {
        teal_error("out of memory");
    }
    return ptr;
}

void *teal_realloc(void *ptr, int num, size_t size)
{
    assert(num >= 0 && size > 0);
    if (!ptr) {
        return teal_alloc(num, size);
    }
    if (num == 0) {
        teal_free(ptr);
        return 0;
    }
    if ((size_t)num > SIZE_MAX / size) {
        teal_error("size overflow");
    }
    void *new = realloc(ptr, (size_t)num * size);
    if (!new) {
        teal_error("out of memory");
    }
    return new;
}

void teal_free(void *ptr)
{
    free(ptr);
}
