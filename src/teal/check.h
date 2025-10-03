#pragma once  // IWYU pragma: always_keep

#define check(expr) ((expr) ? (void)0 : check_fail(__FILE__, __LINE__, __func__, #expr))

void check_fail(const char *file, long line, const char *func, const char *expr);
