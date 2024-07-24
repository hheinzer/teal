#ifndef UTILS_H
#define UTILS_H

#define EPS 1e-8
#define KEYFMT "%30s"

#define ALIAS(a, b) __typeof__(*b) *a = b

#define error(...) x__error(__FILE__, __LINE__, __func__, __VA_ARGS__)
[[noreturn]] void x__error(const char *file, unsigned int line, const char *func,
                           const char *format, ...);

#ifndef NDEBUG
#define ensure(expr) ((expr) ? (void)0 : x__ensure(#expr, __FILE__, __LINE__, __func__))
#else
#define ensure(expr) ((void)0)
#endif
[[noreturn]] void x__ensure(const char *expr, const char *file, unsigned int line,
                            const char *func);

#define sq(a) _Generic(a, long: x__sq_long, double: x__sq_double)(a)
long x__sq_long(long a);
double x__sq_double(double a);

#define min(a, b) _Generic(a, long: x__min_long, double: x__min_double)(a, b)
long x__min_long(long a, long b);
double x__min_double(double a, double b);

#define max(a, b) _Generic(a, long: x__max_long, double: x__max_double)(a, b)
long x__max_long(long a, long b);
double x__max_double(double a, double b);

int isclose(double a, double b);

char sizefmt(double *size);

#endif
