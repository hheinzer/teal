#pragma once

#define EPS 1e-8
#define KEYFMT "%30s"

/* Create an alias 'a' for variable 'b'. */
#define ALIAS(a, b) __typeof__(*b) *a = b

/* Print error message '...' and abort the program. */
#define error(...) x__error(__FILE__, __LINE__, __func__, __VA_ARGS__)
[[gnu::noreturn]] void x__error(const char *file, unsigned int line, const char *func,
                                const char *format, ...);

/* Ensure that expression 'expr' is true, otherwise abort the program. */
#ifndef NDEBUG
#define ensure(expr) ((expr) ? (void)0 : x__ensure(#expr, __FILE__, __LINE__, __func__))
#else
#define ensure(expr) ((void)0)
#endif
[[gnu::noreturn]] void x__ensure(const char *expr, const char *file, unsigned int line,
                                 const char *func);

/* Compute 'a' to the power of 2. */
#define sq(a) _Generic(a, long: x__sq_long, double: x__sq_double)(a)
long x__sq_long(long a);
double x__sq_double(double a);

/* Compute minimum of 'a' and 'b'. */
#define min(a, b) _Generic(a, long: x__min_long, double: x__min_double)(a, b)
long x__min_long(long a, long b);
double x__min_double(double a, double b);

/* Compute maximum of 'a' and 'b'. */
#define max(a, b) _Generic(a, long: x__max_long, double: x__max_double)(a, b)
long x__max_long(long a, long b);
double x__max_double(double a, double b);

/* Determine if 'a' and 'b' are closer than the floating-point tolerance. */
int isclose(double a, double b);

/* Make 'size' human readable (< 1000) and return the corresponding prefix. */
char sizefmt(double *size);
