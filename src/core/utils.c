#include "utils.h"

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "teal.h"

void x__error(const char *file, unsigned int line, const char *func, const char *format, ...)
{
    char msg[128];
    va_list args;
    va_start(args, format);
    vsprintf(msg, format, args);
    va_end(args);
    fprintf(stderr, "[%d] Error: %s (%s:%u: %s)\n", teal.rank, msg, file, line, func);
    int flag;
    MPI_Initialized(&flag);
    if (flag) MPI_Abort(teal.comm, EXIT_FAILURE);
    abort();
}

void x__ensure(const char *expr, const char *file, unsigned int line, const char *func)
{
    fprintf(stderr, "[%d] Assertion '%s' failed (%s:%u: %s)\n", teal.rank, expr, file, line, func);
    int flag;
    MPI_Initialized(&flag);
    if (flag) MPI_Abort(teal.comm, EXIT_FAILURE);
    abort();
}

long x__sq_long(long a)
{
    return a * a;
}
double x__sq_double(double a)
{
    return a * a;
}

long x__min_long(long a, long b)
{
    return a < b ? a : b;
}
double x__min_double(double a, double b)
{
    return a < b ? a : b;
}

long x__max_long(long a, long b)
{
    return a > b ? a : b;
}
double x__max_double(double a, double b)
{
    return a > b ? a : b;
}

int isclose(double a, double b)
{
    return fabs(a - b) < EPS;
}

char sizefmt(double *size)
{
    const char *mod = "\0KMGTPE";  // exabytes might be a bit optimistic...
    while (*size > 1000) {
        *size /= 1000;
        mod++;
    }
    return *mod;
}
