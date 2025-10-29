#include "utils.h"

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "assert.h"
#include "option.h"
#include "sync.h"

void println(const char *fmt, ...)
{
    assert(fmt);
    if (sync.rank == 0 && !option.quiet) {
        va_list args;
        va_start(args, fmt);
        vfprintf(stdout, fmt, args);
        va_end(args);
        fputc('\n', stdout);
        fflush(stdout);
    }
}

void verbose(const char *fmt, ...)
{
    assert(fmt);
    if (sync.rank == 0 && !option.quiet && option.verbose) {
        va_list args;
        va_start(args, fmt);
        vfprintf(stdout, fmt, args);
        va_end(args);
        fputc('\n', stdout);
        fflush(stdout);
    }
}

void error(const char *fmt, ...)
{
    assert(fmt);
    if (sync.rank == 0) {
        fputs("ERROR: ", stderr);
        va_list args;
        va_start(args, fmt);
        vfprintf(stderr, fmt, args);
        va_end(args);
        fputc('\n', stderr);
        fflush(stderr);
    }
    sync_abort();
}

scalar pow2(scalar val)
{
    return val * val;
}

scalar pow3(scalar val)
{
    return val * val * val;
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
    if (lhs == rhs) {
        return true;
    }
    if (!isfinite(lhs) || !isfinite(rhs)) {
        return false;
    }
    static const scalar atol = (sizeof(scalar) == sizeof(float) ? 1e-6 : 1e-12);
    static const scalar rtol = (sizeof(scalar) == sizeof(float) ? 1e-3 : 1e-10);
    return fabs(lhs - rhs) <= fmax(atol, rtol * fmax(fabs(lhs), fabs(rhs)));
}

int cmp_long(const void *lhs_, const void *rhs_)
{
    assert(lhs_ && rhs_);
    const long *lhs = lhs_;
    const long *rhs = rhs_;
    return cmp_asc(*lhs, *rhs);
}

int cmp_scalar(const void *lhs_, const void *rhs_)
{
    assert(lhs_ && rhs_);
    const scalar *lhs = lhs_;
    const scalar *rhs = rhs_;
    return isclose(*lhs, *rhs) ? 0 : cmp_asc(*lhs, *rhs);
}

int cmp_vector(const void *lhs_, const void *rhs_)
{
    assert(lhs_ && rhs_);
    const vector *lhs = lhs_;
    const vector *rhs = rhs_;
    int cmp = cmp_scalar(&lhs->x, &rhs->x);
    if (cmp) {
        return cmp;
    }
    if ((cmp = cmp_scalar(&lhs->y, &rhs->y))) {
        return cmp;
    }
    return cmp_scalar(&lhs->z, &rhs->z);
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

void swap_bytes(void *lhs_, void *rhs_, long size)
{
    assert((lhs_ || rhs_) ? (lhs_ && rhs_ && size >= 0) : size == 0);
    char *lhs = lhs_;
    char *rhs = rhs_;
    for (long i = 0; i < size; i++) {
        char swap = lhs[i];
        lhs[i] = rhs[i];
        rhs[i] = swap;
    }
}

static const char *suffix = "\0KMGTPE";  // ready for exascale computing
static const long base = 1000;

scalar str_to_size(const char *str)
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
    return size;
}

void size_to_str(char *str, scalar size)
{
    assert(str && size >= 0);
    long idx = 0;
    while (size >= base && suffix[idx + 1]) {
        size /= base;
        idx += 1;
    }
    sprintf(str, "%.4g%c", size, suffix[idx]);
}
