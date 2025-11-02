#pragma once

#include <stdint.h>

#include "teal.h"

#define NON_NULL ((void *)(uintptr_t)1)  // NOLINT(performance-no-int-to-ptr)

#define countof(a) (sizeof(a) / sizeof(*(a)))

#define opt(s) (((s) && (s)[0]) ? (s) : "-")

#define cmp_asc(l, r) (((l) > (r)) - ((l) < (r)))
#define cmp_dsc(l, r) cmp_asc(r, l)

void println(const char *fmt, ...) __attribute((format(printf, 1, 2)));
void verbose(const char *fmt, ...) __attribute((format(printf, 1, 2)));
void error(const char *fmt, ...) __attribute((format(printf, 1, 2), noreturn));

scalar pow2(scalar val);
scalar pow3(scalar val);

long lmin(long lhs, long rhs);
long lmax(long lhs, long rhs);

bool isclose(scalar lhs, scalar rhs);

int cmp_long(const void *lhs_, const void *rhs_);
int cmp_scalar(const void *lhs_, const void *rhs_);
int cmp_vector(const void *lhs_, const void *rhs_);

void swap_long(long *lhs, long *rhs);
void swap_scalar(scalar *lhs, scalar *rhs);
void swap_vector(vector *lhs, vector *rhs);
void swap_bytes(void *lhs_, void *rhs_, long size);

scalar str_to_size(const char *str);
void size_to_str(char *str, scalar size);
