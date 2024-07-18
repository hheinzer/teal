#include "array.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "utils.h"

static int longcmp(const void *a, const void *b);
static int doublecmp(const void *a, const void *b);

void x__array_print_long(const long *a, long n)
{
    printf("%ld", a[0]);
    for (long i = 1; i < n; ++i) printf(" %ld", a[i]);
}
void x__array_print_double(const double *a, long n)
{
    printf("%f", a[0]);
    for (long i = 1; i < n; ++i) printf(" %f", a[i]);
}

long x__array_min_long(const long *a, long n)
{
    long min = a[0];
    for (long i = 1; i < n; ++i) min = min(min, a[i]);
    return min;
}
double x__array_min_double(const double *a, long n)
{
    double min = a[0];
    for (long i = 1; i < n; ++i) min = min(min, a[i]);
    return min;
}

long x__array_max_long(const long *a, long n)
{
    long max = a[0];
    for (long i = 1; i < n; ++i) max = max(max, a[i]);
    return max;
}
double x__array_max_double(const double *a, long n)
{
    double max = a[0];
    for (long i = 1; i < n; ++i) max = max(max, a[i]);
    return max;
}

long x__array_sum_long(const long *a, long n)
{
    long sum = a[0];
    for (long i = 1; i < n; ++i) sum += a[i];
    return sum;
}
double x__array_sum_double(const double *a, long n)
{
    double sum = a[0];
    for (long i = 1; i < n; ++i) sum += a[i];
    return sum;
}

long x__array_product_long(const long *a, long n)
{
    long product = a[0];
    for (long i = 1; i < n; ++i) product *= a[i];
    return product;
}
double x__array_product_double(const double *a, long n)
{
    double product = a[0];
    for (long i = 1; i < n; ++i) product *= a[i];
    return product;
}

long x__array_count_long(const long *a, long n, long v)
{
    long count = 0;
    for (long i = 0; i < n; ++i)
        if (a[i] == v) count += 1;
    return count;
}
long x__array_count_double(const double *a, long n, long v)
{
    long count = 0;
    for (long i = 0; i < n; ++i)
        if (isclose(a[i], v)) count += 1;
    return count;
}

long x__array_argmin_long(const long *a, long n)
{
    long arg = 0;
    long min = a[0];
    for (long i = 1; i < n; ++i)
        if (a[i] < min) arg = i;
    return arg;
}
long x__array_argmin_double(const double *a, long n)
{
    long arg = 0;
    double min = a[0];
    for (long i = 1; i < n; ++i)
        if (a[i] < min) arg = i;
    return arg;
}

long x__array_argmax_long(const long *a, long n)
{
    long arg = 0;
    long max = a[0];
    for (long i = 1; i < n; ++i)
        if (a[i] > max) arg = i;
    return arg;
}
long x__array_argmax_double(const double *a, long n)
{
    long arg = 0;
    double max = a[0];
    for (long i = 1; i < n; ++i)
        if (a[i] > max) arg = i;
    return arg;
}

long *x__array_sort_long(long *a, long n)
{
    qsort(a, n, sizeof(*a), longcmp);
    return a;
}
double *x__array_sort_double(double *a, long n)
{
    qsort(a, n, sizeof(*a), doublecmp);
    return a;
}

double array_dot(const double *a, const double *b, long n)
{
    double dot = a[0] * b[0];
    for (long i = 1; i < n; ++i) dot += a[i] * b[i];
    return dot;
}

double array_norm(const double *a, long n)
{
    return sqrt(array_dot(a, a, n));
}

double *array_normalize(double *a, long n)
{
    const double norm = array_norm(a, n);
    for (long i = 0; i < n; ++i) a[i] /= norm;
    return a;
}

double *array_cross(const double *a, const double *b, double *c)
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
    return c;
}

static int longcmp(const void *a, const void *b)
{
    const long *aa = a;
    const long *bb = b;
    if (*aa < *bb) return -1;
    if (*aa > *bb) return 1;
    return 0;
}
static int doublecmp(const void *a, const void *b)
{
    const double *aa = a;
    const double *bb = b;
    if (*aa < *bb) return -1;
    if (*aa > *bb) return 1;
    return 0;
}
