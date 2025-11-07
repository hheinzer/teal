#include "check.h"

#include <stdio.h>

#include "sync.h"

void x__check_fail(const char *file, int line, const char *func, const char *expr)
{
    fprintf(stderr, "[%d] %s:%d: %s: Check `%s` failed.\n", sync.rank, file, line, func, expr);
}
