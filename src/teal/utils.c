#include "utils.h"

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "option.h"
#include "sync.h"

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

int cmp_int(const void *lhs_, const void *rhs_)
{
    assert(lhs_ && rhs_);
    const int *lhs = lhs_;
    const int *rhs = rhs_;
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
    cmp = cmp_scalar(&lhs->y, &rhs->y);
    if (cmp) {
        return cmp;
    }
    return cmp_scalar(&lhs->z, &rhs->z);
}

void swap_int(int *lhs, int *rhs)
{
    assert(lhs && rhs);
    int swap = *lhs;
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

void swap_bytes(void *lhs_, void *rhs_, int size)
{
    assert((lhs_ || rhs_) ? (lhs_ && rhs_ && size >= 0) : (size == 0));
    char *lhs = lhs_;
    char *rhs = rhs_;
    for (int i = 0; i < size; i++) {
        char swap = lhs[i];
        lhs[i] = rhs[i];
        rhs[i] = swap;
    }
}

static const char *suffix = "\0KMGTPE";  // ready for exascale computing
static const int base = 1000;

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
    int idx = 0;
    while (size >= base && suffix[idx + 1]) {
        size /= base;
        idx += 1;
    }
    sprintf(str, "%.4g%c", size, suffix[idx]);
}

void seconds_to_str(char *str, scalar seconds)
{
    assert(str && seconds >= 0);
    static const int seconds_per_minute = 60;
    static const int seconds_per_hour = 60 * seconds_per_minute;
    static const int seconds_per_day = 24 * seconds_per_hour;

    int days = seconds / seconds_per_day;
    seconds -= days * seconds_per_day;

    int hours = seconds / seconds_per_hour;
    seconds -= hours * seconds_per_hour;

    int minutes = seconds / seconds_per_minute;
    seconds -= minutes * seconds_per_minute;

    int len = 0;
    if (days > 0) {
        len += sprintf(str + len, "%dd", days);
    }
    if (hours > 0) {
        len += sprintf(str + len, "%s%dh", len ? " " : "", hours);
    }
    if (minutes > 0) {
        len += sprintf(str + len, "%s%dm", len ? " " : "", minutes);
    }
    sprintf(str + len, "%s%gs", len ? " " : "", seconds);
}
