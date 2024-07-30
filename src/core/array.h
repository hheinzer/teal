#pragma once

/* Print 'n' elements of array 'a'. */
#define array_print(a, n) \
    _Generic(*a, long: x__array_print_long, double: x__array_print_double)(a, n)
void x__array_print_long(const long *a, long n);
void x__array_print_double(const double *a, long n);

/* Compute minimum of 'n' elements of array 'a'. */
#define array_min(a, n) _Generic(*a, long: x__array_min_long, double: x__array_min_double)(a, n)
long x__array_min_long(const long *a, long n);
double x__array_min_double(const double *a, long n);

/* Compute maximum of 'n' elements of array 'a'. */
#define array_max(a, n) _Generic(*a, long: x__array_max_long, double: x__array_max_double)(a, n)
long x__array_max_long(const long *a, long n);
double x__array_max_double(const double *a, long n);

/* Compute sum of 'n' elements of array 'a'. */
#define array_sum(a, n) _Generic(*a, long: x__array_sum_long, double: x__array_sum_double)(a, n)
long x__array_sum_long(const long *a, long n);
double x__array_sum_double(const double *a, long n);

/* Compute product of 'n' elements of array 'a'. */
#define array_product(a, n) \
    _Generic(*a, long: x__array_product_long, double: x__array_product_double)(a, n)
long x__array_product_long(const long *a, long n);
double x__array_product_double(const double *a, long n);

/* Count number of occurrences of 'v' in 'n' elements of array 'a'. */
#define array_count(a, n, v) \
    _Generic(*a, long: x__array_count_long, double: x__array_count_double)(a, n, v)
long x__array_count_long(const long *a, long n, long v);
long x__array_count_double(const double *a, long n, long v);

/* Determine if 'v' is present in 'n' elements of array 'a'. */
#define array_contains(a, n, v) (array_count(a, n, v) > 0)

/* Compute index of minimum in 'n' elements of array 'a'. */
#define array_argmin(a, n) \
    _Generic(*a, long: x__array_argmin_long, double: x__array_argmin_double)(a, n)
long x__array_argmin_long(const long *a, long n);
long x__array_argmin_double(const double *a, long n);

/* Compute index of maximum in 'n' elements of array 'a'. */
#define array_argmax(a, n) \
    _Generic(*a, long: x__array_argmax_long, double: x__array_argmax_double)(a, n)
long x__array_argmax_long(const long *a, long n);
long x__array_argmax_double(const double *a, long n);

/* Sort 'n' elements of array 'a'. */
#define array_sort(a, n) _Generic(*a, long: x__array_sort_long, double: x__array_sort_double)(a, n)
long *x__array_sort_long(long *a, long n);
double *x__array_sort_double(double *a, long n);

/* Compute dot product of 'n' elements of arrays 'a' and 'b'. */
double array_dot(const double *a, const double *b, long n);

/* Compute L2 norm of 'n' elements of array 'a'. */
double array_norm(const double *a, long n);

/* Normalize 'n' elements of array 'a'. */
double *array_normalize(double *a, long n);

/* Compute cross product 'c' of arrays 'a' and 'b'. Each array must be at least 3 elements long. */
double *array_cross(const double *a, const double *b, double *c);
