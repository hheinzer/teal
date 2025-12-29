#include <assert.h>
#include <limits.h>
#include <stdlib.h>

#include "teal.h"

void *teal_realloc(void *ptr, long num, long size)
{
    assert(num >= 0 && size > 0);

    if (!ptr) {
        return teal_alloc(num, size);
    }
    if (num == 0) {
        teal_free(ptr);
        return 0;
    }

    if (num > LONG_MAX / size) {
        teal_error("allocation size overflow");
    }

    void *new = realloc(ptr, num * size);
    if (!new) {
        teal_error("out of memory");
    }

    return new;
}
