#include "array.h"

#include <math.h>
#include <stdio.h>

#include "global.h"
#include "utils.h"

long array_min_long(const long *a, const long n) {
    long min = a[0];
    for (long i = 1; i < n; ++i) min = MIN(min, a[i]);
    return min;
}
double array_min_double(const double *a, const long n) {
    double min = a[0];
    for (long i = 1; i < n; ++i) min = MIN(min, a[i]);
    return min;
}

long array_max_long(const long *a, const long n) {
    long max = a[0];
    for (long i = 1; i < n; ++i) max = MAX(max, a[i]);
    return max;
}
double array_max_double(const double *a, const long n) {
    double max = a[0];
    for (long i = 1; i < n; ++i) max = MAX(max, a[i]);
    return max;
}

long array_sum_long(const long *a, const long n) {
    long sum = 0;
    for (long i = 0; i < n; ++i) sum += a[i];
    return sum;
}
double array_sum_double(const double *a, const long n) {
    double sum = 0;
    for (long i = 0; i < n; ++i) sum += a[i];
    return sum;
}

long array_product_long(const long *a, const long n) {
    long product = 1;
    for (long i = 0; i < n; ++i) product *= a[i];
    return product;
}
double array_product_double(const double *a, const long n) {
    double product = 1;
    for (long i = 0; i < n; ++i) product *= a[i];
    return product;
}

long *array_min_s_long(const long *a, const long n, const long s, long *min) {
    for (long j = 0; j < s; ++j) min[j] = a[0 * s + j];
    for (long i = 1; i < n; ++i)
        for (long j = 0; j < s; ++j) min[j] = MIN(min[j], a[i * s + j]);
    return min;
}
double *array_min_s_double(const double *a, const long n, const long s, double *min) {
    for (long j = 0; j < s; ++j) min[j] = a[0 * s + j];
    for (long i = 1; i < n; ++i)
        for (long j = 0; j < s; ++j) min[j] = MIN(min[j], a[i * s + j]);
    return min;
}

long *array_max_s_long(const long *a, const long n, const long s, long *max) {
    for (long j = 0; j < s; ++j) max[j] = a[0 * s + j];
    for (long i = 1; i < n; ++i)
        for (long j = 0; j < s; ++j) max[j] = MAX(max[j], a[i * s + j]);
    return max;
}
double *array_max_s_double(const double *a, const long n, const long s, double *max) {
    for (long j = 0; j < s; ++j) max[j] = a[0 * s + j];
    for (long i = 1; i < n; ++i)
        for (long j = 0; j < s; ++j) max[j] = MAX(max[j], a[i * s + j]);
    return max;
}

long array_argmin_long(const long *a, const long n) {
    long argmin = 0;
    long min = a[0];
    for (long i = 1; i < n; ++i)
        if (a[i] < min) argmin = i;
    return argmin;
}
long array_argmin_double(const double *a, const long n) {
    long argmin = 0;
    double min = a[0];
    for (long i = 1; i < n; ++i)
        if (a[i] < min) argmin = i;
    return argmin;
}

long array_argmax_long(const long *a, const long n) {
    long argmax = 0;
    long max = a[0];
    for (long i = 1; i < n; ++i)
        if (a[i] > max) argmax = i;
    return argmax;
}
long array_argmax_double(const double *a, const long n) {
    long argmax = 0;
    double max = a[0];
    for (long i = 1; i < n; ++i)
        if (a[i] > max) argmax = i;
    return argmax;
}

bool array_contains_long(const long *a, const long n, const long v) {
    for (long i = 0; i < n; ++i)
        if (a[i] == v) return true;
    return false;
}
bool array_contains_double(const double *a, const long n, const double v) {
    for (long i = 0; i < n; ++i)
        if (fabs(a[i] - v) < EPS) return true;
    return false;
}

long array_count_long(const long *a, const long n, const long v) {
    long count = 0;
    for (long i = 0; i < n; ++i)
        if (a[i] == v) count += 1;
    return count;
}
long array_count_double(const double *a, const long n, const double v) {
    long count = 0;
    for (long i = 0; i < n; ++i)
        if (fabs(a[i] - v) < EPS) count += 1;
    return count;
}

double array_norm2(const double *a, const long n) {
    double norm2 = 0;
    for (long i = 0; i < n; ++i) norm2 += a[i] * a[i];
    return sqrt(norm2);
}

double *array_norm2_s(const double *a, const long n, const long s, double *norm2) {
    for (long j = 0; j < s; ++j) norm2[j] = 0;
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < s; ++j) norm2[j] += a[i * s + j] * a[i * s + j];
    for (long j = 0; j < s; ++j) norm2[j] = sqrt(norm2[j]);
    return norm2;
}

void array_normalize(double *a, const long n) {
    double norm2 = 0;
    for (long i = 0; i < n; ++i) norm2 += a[i] * a[i];
    norm2 = sqrt(norm2);
    for (long i = 0; i < n; ++i) a[i] /= norm2;
}

double array_dot(const double *a, const double *b, const long n) {
    double dot = 0;
    for (long i = 0; i < n; ++i) dot += a[i] * b[i];
    return dot;
}

void array_print_long(const char *pre, const long *a, const long n, const char *post) {
    if (pre) printf("%s", pre);
    if (n) printf("%ld", a[0]);
    for (long i = 1; i < n; ++i) printf(" %ld", a[i]);
    if (post) printf("%s", post);
}
void array_print_double(const char *pre, const double *a, const long n, const char *post) {
    if (pre) printf("%s", pre);
    if (n) printf("%g", a[0]);
    for (long i = 1; i < n; ++i) printf(" %g", a[i]);
    if (post) printf("%s", post);
}

void array_fprint_long(FILE *file, const long *a, const long n) {
    for (long i = 0; i < n; ++i) fprintf(file, "%ld ", a[i]);
    fputc('\n', file);
}
void array_fprint_double(FILE *file, const double *a, const long n) {
    for (long i = 0; i < n; ++i) fprintf(file, "%g ", a[i]);
    fputc('\n', file);
}
