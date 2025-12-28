#include <assert.h>
#include <stdarg.h>
#include <stdio.h>

#include "sync.h"
#include "teal.h"

void teal_print(const char *fmt, ...)
{
    assert(fmt);
    if (sync.rank == 0 && !teal.quiet) {
        char msg[128];

        va_list args;
        va_start(args, fmt);
        int len = vsprintf(msg, fmt, args);
        va_end(args);

        msg[len] = '\n';
        msg[len + 1] = 0;
        fputs(msg, stdout);
        fflush(stdout);
    }
}
