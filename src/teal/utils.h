#pragma once

#include <stdint.h>

#include "teal.h"

#define NON_NULL ((void *)(uintptr_t)1)  // NOLINT(performance-no-int-to-ptr)

#define countof(a) (sizeof(a) / sizeof(*(a)))

#define cmp_asc(l, r) (((l) > (r)) - ((l) < (r)))
#define cmp_dsc(l, r) cmp_asc(r, l)

void println(const char *fmt, ...) __attribute((format(printf, 1, 2)));
void verbose(const char *fmt, ...) __attribute((format(printf, 1, 2)));
void error(const char *fmt, ...) __attribute((format(printf, 1, 2), noreturn));

scalar sq(scalar val);

long lmin(long lhs, long rhs);
long lmax(long lhs, long rhs);

bool isclose(scalar lhs, scalar rhs);

int cmp_long(const void *lhs_, const void *rhs_);
int cmp_scalar(const void *lhs_, const void *rhs_);
int cmp_vector(const void *lhs_, const void *rhs_);

void lswap(long *lhs, long *rhs);
void fswap(scalar *lhs, scalar *rhs);
void vswap(vector *lhs, vector *rhs);
void swap_bytes(void *lhs_, void *rhs_, long size);

scalar str_to_size(const char *str);
void size_to_str(char *str, scalar size);
