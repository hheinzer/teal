#include <assert.h>
#include <stdint.h>
#include <string.h>

#include "teal.h"
#include "utils.h"

void *memdup(const void *ptr, int num, size_t size)
{
    assert((ptr || num == 0) && num >= 0 && size > 0);
    if (num == 0) {
        return 0;
    }
    void *dup = teal_alloc(num, size);
    assert((size_t)num <= SIZE_MAX / size);
    memcpy(dup, ptr, (size_t)num * size);
    return dup;
}
