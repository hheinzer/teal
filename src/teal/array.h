#pragma once

#define array_print(a, n) \
    _Generic(*a, long: x__array_print_long, double: x__array_print_double)(a, n)
void x__array_print_long(const long *a, long n);
void x__array_print_double(const double *a, long n);

#define array_min(a, n) _Generic(*a, long: x__array_min_long, double: x__array_min_double)(a, n)
long x__array_min_long(const long *a, long n);
double x__array_min_double(const double *a, long n);

#define array_max(a, n) _Generic(*a, long: x__array_max_long, double: x__array_max_double)(a, n)
long x__array_max_long(const long *a, long n);
double x__array_max_double(const double *a, long n);

#define array_sum(a, n) _Generic(*a, long: x__array_sum_long, double: x__array_sum_double)(a, n)
long x__array_sum_long(const long *a, long n);
double x__array_sum_double(const double *a, long n);

#define array_product(a, n) \
    _Generic(*a, long: x__array_product_long, double: x__array_product_double)(a, n)
long x__array_product_long(const long *a, long n);
double x__array_product_double(const double *a, long n);

#define array_count(a, n, v) \
    _Generic(*a, long: x__array_count_long, double: x__array_count_double)(a, n, v)
long x__array_count_long(const long *a, long n, long v);
long x__array_count_double(const double *a, long n, double v);

#define array_contains(a, n, v) \
    _Generic(*a, long: x__array_contains_long, double: x__array_contains_double)(a, n, v)
int x__array_contains_long(const long *a, long n, long v);
int x__array_contains_double(const double *a, long n, double v);

#define array_argmin(a, n) \
    _Generic(*a, long: x__array_argmin_long, double: x__array_argmin_double)(a, n)
long x__array_argmin_long(const long *a, long n);
long x__array_argmin_double(const double *a, long n);

#define array_argmax(a, n) \
    _Generic(*a, long: x__array_argmax_long, double: x__array_argmax_double)(a, n)
long x__array_argmax_long(const long *a, long n);
long x__array_argmax_double(const double *a, long n);

#define array_sort(a, n) _Generic(*a, long: x__array_sort_long, double: x__array_sort_double)(a, n)
long *x__array_sort_long(long *a, long n);
double *x__array_sort_double(double *a, long n);

double array_dot(const double *a, const double *b, long n);

double array_norm(const double *a, long n);

double *array_normalize(double *a, long n);

double *array_cross(const double a[3], const double b[3], double c[3]);
