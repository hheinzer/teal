#pragma once

#include "teal.h"

#define NON_NULL ((void *)(number)1)  // NOLINT(performance-no-int-to-ptr)

#define countof(arr) (sizeof(arr) / sizeof(*(arr)))

#define optional(str) (((str) && *(str)) ? (str) : "none")

#define cmp_asc(lhs, rhs) (((lhs) > (rhs)) - ((lhs) < (rhs)))
#define cmp_dsc(lhs, rhs) cmp_asc(rhs, lhs)

void println(const char *fmt, ...) __attribute((format(printf, 1, 2)));
void verbose(const char *fmt, ...) __attribute((format(printf, 1, 2)));
void error(const char *fmt, ...) __attribute((format(printf, 1, 2), noreturn));

scalar pow2(scalar val);
scalar pow3(scalar val);

number lmin(number lhs, number rhs);
number lmax(number lhs, number rhs);

bool isclose(scalar lhs, scalar rhs);

int cmp_number(const void *lhs_, const void *rhs_);
int cmp_scalar(const void *lhs_, const void *rhs_);
int cmp_vector(const void *lhs_, const void *rhs_);

void swap_number(number *lhs, number *rhs);
void swap_scalar(scalar *lhs, scalar *rhs);
void swap_vector(vector *lhs, vector *rhs);
void swap_bytes(void *lhs_, void *rhs_, number size);

scalar str_to_size(const char *str);
void size_to_str(char *str, scalar size);
