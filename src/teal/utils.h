#pragma once

#include <stdint.h>

#include "teal.h"

#define unused(val) ((void)(val))

#define countof(arr) (sizeof(arr) / sizeof(*(arr)))

#define cmp_asc(lhs, rhs) (((lhs) > (rhs)) - ((lhs) < (rhs)))
#define cmp_dsc(lhs, rhs) cmp_asc(rhs, lhs)

scalar pow2(scalar val);
scalar pow3(scalar val);

int lmin(int lhs, int rhs);
int lmax(int lhs, int rhs);

bool isclose(scalar lhs, scalar rhs);

void println(const char *fmt, ...) __attribute((format(printf, 1, 2)));
void verbose(const char *fmt, ...) __attribute((format(printf, 1, 2)));
void error(const char *fmt, ...) __attribute((format(printf, 1, 2), noreturn));

int cmp_int(const void *lhs_, const void *rhs_);
int cmp_scalar(const void *lhs_, const void *rhs_);
int cmp_vector(const void *lhs_, const void *rhs_);

void swap_int(int *lhs, int *rhs);
void swap_scalar(scalar *lhs, scalar *rhs);
void swap_vector(vector *lhs, vector *rhs);
void swap_bytes(void *lhs_, void *rhs_, int size);

scalar str_to_size(const char *str);
void size_to_str(char *str, scalar size);

void seconds_to_str(char *str, scalar seconds);
