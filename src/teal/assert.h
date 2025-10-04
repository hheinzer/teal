#pragma once  // IWYU pragma: always_keep

#define STATIC_ASSERT(expr) extern int STATIC_ASSERT[(expr) ? 1 : -1]

#ifndef __has_builtin
#define __has_builtin(x) 0
#endif

#if __has_builtin(__builtin_expect)
#define expect(expr, cond) __builtin_expect(expr, cond)
#else
#define expect(expr, cond) (expr)
#endif

#if __has_builtin(__builtin_unreachable)
#define unreachable() __builtin_unreachable()
#else
__attribute((noreturn)) static inline void unreachable(void)
{
}
#endif

#ifdef assert
#undef assert
#endif

#ifdef NDEBUG
#define assert(expr) (expect(!!(expr), 1) ? (void)0 : unreachable())
#else
#define assert(expr) ((expr) ? (void)0 : x__assert_fail(__FILE__, __LINE__, __func__, #expr))
#endif

void x__assert_fail(const char *file, long line, const char *func, const char *expr)
    __attribute((noreturn));
