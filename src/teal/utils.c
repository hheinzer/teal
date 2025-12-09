#include "utils.h"

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "option.h"
#include "sync.h"

scalar sq(scalar val)
{
    return val * val;
}

scalar cb(scalar val)
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

bool is_close(scalar lhs, scalar rhs)
{
    assert(isfinite(lhs) && isfinite(rhs));
    static const scalar atol = (sizeof(scalar) == sizeof(float) ? 1e-6 : 1e-8);
    static const scalar rtol = (sizeof(scalar) == sizeof(float) ? 1e-3 : 1e-5);
    return fabs(lhs - rhs) <= fmax(atol, rtol * fmax(fabs(lhs), fabs(rhs)));
}

bool is_less(scalar lhs, scalar rhs)
{
    return !is_close(lhs, rhs) && lhs < rhs;
}

bool is_greater(scalar lhs, scalar rhs)
{
    return !is_close(lhs, rhs) && lhs > rhs;
}

bool is_close_or_less(scalar lhs, scalar rhs)
{
    return is_close(lhs, rhs) || lhs < rhs;
}

bool is_close_or_greater(scalar lhs, scalar rhs)
{
    return is_close(lhs, rhs) || lhs > rhs;
}

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
    MPI_Abort(sync.comm, EXIT_FAILURE);
    abort();
}

int cmp_long(const void *lhs_, const void *rhs_)
{
    assert(lhs_ && rhs_);
    const long *lhs = lhs_;
    const long *rhs = rhs_;
    return (*lhs > *rhs) - (*lhs < *rhs);
}

int cmp_scalar(const void *lhs_, const void *rhs_)
{
    assert(lhs_ && rhs_);
    const scalar *lhs = lhs_;
    const scalar *rhs = rhs_;
    return is_close(*lhs, *rhs) ? 0 : (*lhs > *rhs) - (*lhs < *rhs);
}

int cmp_vector(const void *lhs_, const void *rhs_)
{
    assert(lhs_ && rhs_);
    const vector *lhs = lhs_;
    const vector *rhs = rhs_;
    long cmp = cmp_scalar(&lhs->x, &rhs->x);
    if (cmp) {
        return cmp;
    }
    cmp = cmp_scalar(&lhs->y, &rhs->y);
    if (cmp) {
        return cmp;
    }
    return cmp_scalar(&lhs->z, &rhs->z);
}

void swap_long(long *lhs, long *rhs)
{
    assert(lhs && rhs);
    long swap = *lhs;
    *lhs = *rhs;
    *rhs = swap;
}

void swap_scalar(scalar *lhs, scalar *rhs)
{
    assert(lhs && rhs);
    scalar swap = *lhs;
    *lhs = *rhs;
    *rhs = swap;
}

void swap_vector(vector *lhs, vector *rhs)
{
    assert(lhs && rhs);
    vector swap = *lhs;
    *lhs = *rhs;
    *rhs = swap;
}

void swap_bytes(void *lhs_, void *rhs_, long size)
{
    assert((lhs_ || rhs_) ? (lhs_ && rhs_ && size >= 0) : (size == 0));
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

void seconds_to_str(char *str, scalar seconds)
{
    assert(str && seconds >= 0);
    enum { SECONDS_PER_MINUTE = 60 };
    enum { SECONDS_PER_HOUR = 60 * SECONDS_PER_MINUTE };
    enum { SECONDS_PER_DAY = 24 * SECONDS_PER_HOUR };

    long days = seconds / SECONDS_PER_DAY;
    seconds -= days * SECONDS_PER_DAY;

    long hours = seconds / SECONDS_PER_HOUR;
    seconds -= hours * SECONDS_PER_HOUR;

    long minutes = seconds / SECONDS_PER_MINUTE;
    seconds -= minutes * SECONDS_PER_MINUTE;

    long len = 0;
    if (days > 0) {
        len += sprintf(str + len, "%ldd", days);
    }
    if (hours > 0) {
        len += sprintf(str + len, "%s%ldh", len ? " " : "", hours);
    }
    if (minutes > 0) {
        len += sprintf(str + len, "%s%ldm", len ? " " : "", minutes);
    }
    sprintf(str + len, "%s%gs", len ? " " : "", seconds);
}
