#pragma once  // IWYU pragma: always_keep

#define STATIC_ASSERT(expr) extern int STATIC_ASSERT[(expr) ? 1 : -1]

#ifdef NDEBUG
#define assert(expr) (__builtin_expect(!!(expr), 1) ? (void)0 : __builtin_unreachable())
#else
#define assert(expr) ((expr) ? (void)0 : x__assert_fail(__FILE__, __LINE__, __func__, #expr))
#endif

void x__assert_fail(const char *file, long line, const char *func, const char *expr)
    __attribute((noreturn));
