#pragma once

#include <stdint.h>

#include "teal.h"

static const int NON_NULL;
#define NON_NULL ((void *)&NON_NULL)

#define countof(a) (sizeof(a) / sizeof(*(a)))

#ifdef NDEBUG
#define assert(expr) (__builtin_expect(!!(expr), 1) ? (void)0 : __builtin_unreachable())
#else
#define assert(expr) ((expr) ? (void)0 : assert_fail(__FILE__, __LINE__, __func__, #expr))
#endif

void assert_fail(const char *file, long line, const char *func, const char *expr)
    __attribute((noreturn));

void print(const char *format, ...) __attribute((format(printf, 1, 2)));

double sq(double val);

long lmin(long lhs, long rhs);
long lmax(long lhs, long rhs);

bool isclose(double lhs, double rhs);

int lcmp(const void *lhs, const void *rhs);
int fcmp(const void *lhs, const void *rhs);

void lswap(long *lhs, long *rhs);
void fswap(double *lhs, double *rhs);
void vswap(vector *lhs, vector *rhs);
void memswap(void *lhs, void *rhs, long size);

uint64_t fnv1a(const void *ptr, long size);

bool fexists(const char *fname);

int strrot(char *dst, const char *src, char sep);

long str2size(const char *str);
void size2str(char *str, double size);
