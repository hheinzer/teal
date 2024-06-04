#ifndef ARRAY_H
#define ARRAY_H

#include <stdio.h>

long array_min_long(const long *a, const long n);
double array_min_double(const double *a, const long n);
#define array_min(a, n)               \
    _Generic(a,                       \
        long *: array_min_long,       \
        double *: array_min_double,   \
        const long *: array_min_long, \
        const double *: array_min_double)(a, n)

long array_max_long(const long *a, const long n);
double array_max_double(const double *a, const long n);
#define array_max(a, n)               \
    _Generic(a,                       \
        long *: array_max_long,       \
        double *: array_max_double,   \
        const long *: array_max_long, \
        const double *: array_max_double)(a, n)

long array_sum_long(const long *a, const long n);
double array_sum_double(const double *a, const long n);
#define array_sum(a, n)               \
    _Generic(a,                       \
        long *: array_sum_long,       \
        double *: array_sum_double,   \
        const long *: array_sum_long, \
        const double *: array_sum_double)(a, n)

long array_product_long(const long *a, const long n);
double array_product_double(const double *a, const long n);
#define array_product(a, n)               \
    _Generic(a,                           \
        long *: array_product_long,       \
        double *: array_product_double,   \
        const long *: array_product_long, \
        const double *: array_product_double)(a, n)

long *array_min_s_long(const long *a, const long n, const long s, long *min);
double *array_min_s_double(const double *a, const long n, const long s, double *min);
#define array_min_s(a, n, s, min)       \
    _Generic(a,                         \
        long *: array_min_s_long,       \
        double *: array_min_s_double,   \
        const long *: array_min_s_long, \
        const double *: array_min_s_double)(a, n, s, min)

long *array_max_s_long(const long *a, const long n, const long s, long *max);
double *array_max_s_double(const double *a, const long n, const long s, double *max);
#define array_max_s(a, n, s, max)       \
    _Generic(a,                         \
        long *: array_max_s_long,       \
        double *: array_max_s_double,   \
        const long *: array_max_s_long, \
        const double *: array_max_s_double)(a, n, s, max)

long array_argmin_long(const long *a, const long n);
long array_argmin_double(const double *a, const long n);
#define array_argmin(a, n)               \
    _Generic(a,                          \
        long *: array_argmin_long,       \
        double *: array_argmin_double,   \
        const long *: array_argmin_long, \
        const double *: array_argmin_double)(a, n)

long array_argmax_long(const long *a, const long n);
long array_argmax_double(const double *a, const long n);
#define array_argmax(a, n)               \
    _Generic(a,                          \
        long *: array_argmax_long,       \
        double *: array_argmax_double,   \
        const long *: array_argmax_long, \
        const double *: array_argmax_double)(a, n)

bool array_contains_long(const long *a, const long n, const long v);
bool array_contains_double(const double *a, const long n, const double v);
#define array_contains(a, n, v)            \
    _Generic(a,                            \
        long *: array_contains_long,       \
        double *: array_contains_double,   \
        const long *: array_contains_long, \
        const double *: array_contains_double)(a, n, v)

long array_count_long(const long *a, const long n, const long v);
long array_count_double(const double *a, const long n, const double v);
#define array_count(a, n, v)            \
    _Generic(a,                         \
        long *: array_count_long,       \
        double *: array_count_double,   \
        const long *: array_count_long, \
        const double *: array_count_double)(a, n, v)

double array_norm2(const double *a, const long n);

double *array_norm2_s(const double *a, const long n, const long s, double *norm2);

void array_normalize(double *a, const long n);

double array_dot(const double *a, const double *b, const long n);

void array_print_long(const long *a, const long n, const char *post);
void array_print_double(const double *a, const long n, const char *post);
#define array_print(a, n, post)         \
    _Generic(a,                         \
        long *: array_print_long,       \
        double *: array_print_double,   \
        const long *: array_print_long, \
        const double *: array_print_double)(a, n, post)

void array_fprint_long(FILE *file, const long *a, const long n);
void array_fprint_double(FILE *file, const double *a, const long n);
#define array_fprint(file, a, n)         \
    _Generic(a,                          \
        long *: array_fprint_long,       \
        double *: array_fprint_double,   \
        const long *: array_fprint_long, \
        const double *: array_fprint_double)(file, a, n)

#endif
