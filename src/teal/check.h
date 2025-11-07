#pragma once  // IWYU pragma: always_keep

#define check(expr) ((expr) ? (void)0 : x__check_fail(__FILE__, __LINE__, __func__, #expr))

void x__check_fail(const char *file, int line, const char *func, const char *expr);
