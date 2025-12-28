#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "sync.h"
#include "teal.h"

void teal_error(const char *fmt, ...)
{
    assert(fmt);

    char msg[128];
    int len = sprintf(msg, "[%d] ERROR: ", sync.rank);

    va_list args;
    va_start(args, fmt);
    len += vsprintf(msg + len, fmt, args);
    va_end(args);

    msg[len] = '\n';
    msg[len + 1] = 0;
    fputs(msg, stderr);
    fflush(stderr);

    teal_exit(EXIT_FAILURE);
}
