#pragma once

#include "teal.h"

#define NON_NULL ((void *)(int)1)  // NOLINT(performance-no-int-to-ptr)

#define count_of(arr) (sizeof(arr) / sizeof(*(arr)))

void println(const char *fmt, ...) __attribute((format(printf, 1, 2)));
void verbose(const char *fmt, ...) __attribute((format(printf, 1, 2)));
void error(const char *fmt, ...) __attribute((format(printf, 1, 2), noreturn));

scalar sq(scalar val);
scalar cb(scalar val);

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

void seconds_to_str(char *str, scalar seconds);
