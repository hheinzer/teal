#include <assert.h>
#include <limits.h>
#include <stdlib.h>

#include "teal.h"

void *teal_alloc(long num, long size)
{
    assert(num >= 0 && size > 0);

    if (num == 0) {
        return 0;
    }

    if (num > LONG_MAX / size) {
        teal_error("allocation size overflow");
    }

    void *ptr = calloc((size_t)num, (size_t)size);
    if (!ptr) {
        teal_error("out of memory");
    }

    return ptr;
}
