#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "teal.h"
#include "teal/sync.h"

void teal_exit(const char *format, ...)
{
    if (sync.rank == 0) {
        va_list args;
        va_start(args, format);
        printf("\n");
        vprintf(format, args);
        printf("\n\n");
        va_end(args);
    }
    teal_finalize();
    exit(EXIT_SUCCESS);
}
