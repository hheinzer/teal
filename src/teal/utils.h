#pragma once

#include <stdint.h>

#include "teal.h"

#define NON_NULL ((void *)(uintptr_t)1)  // NOLINT(performance-no-int-to-ptr)

#define countof(a) (sizeof(a) / sizeof(*(a)))

void print(const char *fmt, ...) __attribute((format(printf, 1, 2)));

scalar sq(scalar val);

long lmin(long lhs, long rhs);
long lmax(long lhs, long rhs);

bool isclose(scalar lhs, scalar rhs);

int lcmp(const void *lhs, const void *rhs);
int fcmp(const void *lhs, const void *rhs);
int vcmp(const void *lhs, const void *rhs);

void lswap(long *lhs, long *rhs);
void fswap(scalar *lhs, scalar *rhs);
void vswap(vector *lhs, vector *rhs);
void memswap(void *lhs, void *rhs, long size);

bool fexists(const char *fname);

scalar str2size(const char *str);
void size2str(char *str, scalar size);
