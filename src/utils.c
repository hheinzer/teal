#include "utils.h"

#include <assert.h>
#include <string.h>

#include "memory.h"

void *utils_memdup(const void *src, const long nmemb, const long size) {
    if (!src) return 0;
    void *dst = memory_calloc(nmemb, size);
    dst = memcpy(dst, src, nmemb * size);
    assert(dst && "memcpy failure");
    return dst;
}

char *utils_strdup(const char *src) {
    return utils_memdup(src, strlen(src) + 1, sizeof(*src));
}

const char *utils_extension(const char *fname) {
    const char *dot = strrchr(fname, '.');
    assert(dot && "file name has no extension");
    return dot + 1;
}

const char *utils_basename(const char *fname) {
    const char *slash = strrchr(fname, '/');
    if (!slash) return fname;
    return slash + 1;
}
