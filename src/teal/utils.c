#include "utils.h"

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "option.h"
#include "sync.h"

void assert_fail(const char *file, long line, const char *func, const char *expr)
{
    fprintf(stderr, "[%d] %s:%ld: %s: Assertion `%s` failed.\n", sync.rank, file, line, func, expr);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    abort();
}

void print(const char *format, ...)
{
    assert(format);
    if (sync.rank == 0 && !option.quiet) {
        va_list args;
        va_start(args, format);
        vprintf(format, args);
        va_end(args);
    }
}

double sq(double val)
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

bool isclose(double lhs, double rhs)
{
    static const double atol = 1.0e-12;
    static const double rtol = 1.0e-10;
    return fabs(lhs - rhs) <= fmax(atol, rtol * fmax(fabs(lhs), fabs(rhs)));
}

int lcmp(const void *lhs, const void *rhs)
{
    assert(lhs && rhs);
    long _lhs = *(long *)lhs;
    long _rhs = *(long *)rhs;
    return (_lhs > _rhs) - (_lhs < _rhs);
}

int fcmp(const void *lhs, const void *rhs)
{
    assert(lhs && rhs);
    double _lhs = *(double *)lhs;
    double _rhs = *(double *)rhs;
    return isclose(_lhs, _rhs) ? 0 : (_lhs > _rhs) - (_lhs < _rhs);
}

void lswap(long *lhs, long *rhs)
{
    assert(lhs && rhs);
    long swap = *lhs;
    *lhs = *rhs;
    *rhs = swap;
}

void fswap(double *lhs, double *rhs)
{
    assert(lhs && rhs);
    double swap = *lhs;
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

uint64_t fnv1a(const void *ptr, long size)
{
    assert(ptr ? size >= 0 : size == 0);
    static const uint64_t basis = 0xcbf29ce484222325;
    static const uint64_t prime = 0x00000100000001b3;
    const uint8_t *byte = ptr;
    uint64_t hash = basis;
    for (long i = 0; i < size; i++) {
        hash ^= byte[i];
        hash *= prime;
    }
    return hash;
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

int strrot(char *dst, const char *src, char sep)
{
    assert(dst && src);
    char *pos = strchr(src, sep);
    return pos ? sprintf(dst, "%s%c%.*s", pos + 1, sep, (int)(pos - src), src) : 0;
}

long str2size(const char *str)
{
    assert(str);
    static const double base = 1000;
    char *end;
    double size = strtod(str, &end);
    switch (*end) {
        case 'K': size *= base; break;
        case 'M': size *= base * base; break;
        case 'G': size *= base * base * base; break;
        case 'T': size *= base * base * base * base; break;
        default: break;
    }
    return ceil(size);
}

void size2str(char *str, double size)
{
    assert(str && size >= 0);
    static const double base = 1000;
    char *prefix = "\0KMGT";
    while (size >= base && prefix[1]) {
        size /= base;
        prefix += 1;
    }
    sprintf(str, "%.4g%c", size, *prefix);
}
