#include "utils.h"

#define UNW_LOCAL_ONLY
#include <libunwind.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "option.h"
#include "sync.h"

enum { BUFLEN = 4096 };

static int append(char *str, long *pos, const char *format, ...)
{
    if (!str || !pos || *pos < 0 || !format) {
        abort();
    }
    if (*pos >= BUFLEN) {
        str[BUFLEN - 1] = 0;
        return 1;
    }

    va_list args;
    va_start(args, format);
    long rem = BUFLEN - *pos;
    long num = vsnprintf(&str[*pos], rem, format, args);
    va_end(args);

    if (num < 0) {
        return -1;
    }
    if (num >= rem) {
        *pos = BUFLEN - 1;
        str[*pos] = 0;
        return 1;
    }
    *pos += num;
    return 0;
}

void assert_fail(const char *file, long line, const char *func, const char *expr)
{
    if (!file || line < 0 || !func || !expr) {
        abort();
    }

    char buf[BUFLEN];
    long pos = 0;

    if (append(buf, &pos, "[%d] %s:%ld: %s: Assertion `%s` failed.\n", sync.rank, file, line, func,
               expr)) {
        goto out;
    }

    unw_context_t ctx;
    unw_cursor_t cur;
    if (!unw_getcontext(&ctx) && !unw_init_local(&cur, &ctx)) {
        long frame = 0;
        while (unw_step(&cur) > 0) {
            string name = "";
            unw_word_t offset = 0;
            if (unw_get_proc_name(&cur, name, sizeof(name), &offset)) {
                strcpy(name, "???");
                offset = 0;
            }
            unw_word_t iptr = 0;
            unw_get_reg(&cur, UNW_REG_IP, &iptr);
            if (append(buf, &pos, "\t %2ld. %-30s (+0x%lx) [0x%lx]\n", frame++, name, offset,
                       iptr)) {
                break;
            }
        }
    }

out:
    if (pos > 0) {
        fwrite(buf, sizeof(*buf), pos, stderr);
        fflush(stderr);
    }

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
    const long *_lhs = lhs;
    const long *_rhs = rhs;
    return (*_lhs > *_rhs) - (*_lhs < *_rhs);
}

int fcmp(const void *lhs, const void *rhs)
{
    assert(lhs && rhs);
    const double *_lhs = lhs;
    const double *_rhs = rhs;
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
    cmp = fcmp(&_lhs->y, &_rhs->y);
    if (cmp) {
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
