#include "utils.h"

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "assert.h"
#include "option.h"
#include "sync.h"

void print(const char *fmt, ...)
{
    assert(fmt);
    if (sync.rank == 0 && !option.quiet) {
        va_list args;
        va_start(args, fmt);
        vprintf(fmt, args);
        va_end(args);
    }
}

scalar sq(scalar val)
{
    return val * val;
}

long lmin(long lhs, long rhs)
{
    return (lhs < rhs) ? lhs : rhs;
}

long lmax(long lhs, long rhs)
{
    return (lhs > rhs) ? lhs : rhs;
}

bool isclose(scalar lhs, scalar rhs)
{
    static const scalar atol = (sizeof(scalar) == sizeof(float) ? 1e-6 : 1e-12);
    static const scalar rtol = (sizeof(scalar) == sizeof(float) ? 1e-3 : 1e-10);
    return fabs(lhs - rhs) <= fmax(atol, rtol * fmax(fabs(lhs), fabs(rhs)));
}

int lcmp(const void *lhs, const void *rhs)
{
    assert(lhs && rhs);
    const long *_lhs = lhs;
    const long *_rhs = rhs;
    return (*_lhs > *_rhs) - (*_lhs < *_rhs);
}

int fcmp(const void *lhs, const void *rhs)
{
    assert(lhs && rhs);
    const scalar *_lhs = lhs;
    const scalar *_rhs = rhs;
    return isclose(*_lhs, *_rhs) ? 0 : (*_lhs > *_rhs) - (*_lhs < *_rhs);
}

int vcmp(const void *lhs, const void *rhs)
{
    assert(lhs && rhs);
    const vector *_lhs = lhs;
    const vector *_rhs = rhs;
    int cmp = fcmp(&_lhs->x, &_rhs->x);
    if (cmp) {
        return cmp;
    }
    if ((cmp = fcmp(&_lhs->y, &_rhs->y))) {
        return cmp;
    }
    return fcmp(&_lhs->z, &_rhs->z);
}

void lswap(long *lhs, long *rhs)
{
    assert(lhs && rhs);
    long swap = *lhs;
    *lhs = *rhs;
    *rhs = swap;
}

void fswap(scalar *lhs, scalar *rhs)
{
    assert(lhs && rhs);
    scalar swap = *lhs;
    *lhs = *rhs;
    *rhs = swap;
}

void vswap(vector *lhs, vector *rhs)
{
    assert(lhs && rhs);
    vector swap = *lhs;
    *lhs = *rhs;
    *rhs = swap;
}

void memswap(void *lhs, void *rhs, long size)
{
    assert((lhs || rhs) ? (lhs && rhs && size >= 0) : size == 0);
    char *_lhs = lhs;
    char *_rhs = rhs;
    for (long i = 0; i < size; i++) {
        char swap = _lhs[i];
        _lhs[i] = _rhs[i];
        _rhs[i] = swap;
    }
}

bool fexists(const char *fname)
{
    assert(fname);
    FILE *file = fopen(fname, "r");
    if (!file) {
        return false;
    }
    fclose(file);
    return true;
}

static const char suffix[] = "\0KMGTPE";  // ready for exascale computing
static const long base = 1000;

long str2size(const char *str)
{
    assert(str);
    char *end;
    scalar size = strtod(str, &end);
    if (*end) {
        char *pos = strchr(&suffix[1], *end);
        if (pos) {
            size *= pow(base, pos - suffix);
        }
    }
    return ceil(size);
}

strbuf size2str(scalar size)
{
    assert(size >= 0);
    long idx = 0;
    while (size >= base && suffix[idx + 1]) {
        size /= base;
        idx += 1;
    }
    strbuf str;
    sprintf(str.buf, "%.4g%c", size, suffix[idx]);
    return str;
}
