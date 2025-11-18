#pragma once

#include <stdint.h>

#include "teal.h"

#define unused(val) ((void)(val))

#define countof(arr) (sizeof(arr) / sizeof(*(arr)))

#define sq(val) ((val) * (val))
#define cb(val) ((val) * (val) * (val))

#define min(lhs, rhs) (((lhs) < (rhs)) ? (lhs) : (rhs))
#define max(lhs, rhs) (((lhs) > (rhs)) ? (lhs) : (rhs))

#define cmp_asc(lhs, rhs) (((lhs) > (rhs)) - ((lhs) < (rhs)))
#define cmp_dsc(lhs, rhs) cmp_asc(rhs, lhs)

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
