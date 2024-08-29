#include "print.h"

#include <stdarg.h>
#include <stdio.h>

void print_key(const char *key, const char *format, ...)
{
    va_list args;
    va_start(args, format);
    printf(" | %30s: ", key);
    vprintf(format, args);
    printf("\n");
    va_end(args);
}

void print_size(double size_bytes)
{
    const char *prefix = "\0KMGTPE";  // exabytes might be a bit optimistic...
    while (size_bytes > 1000) {
        size_bytes /= 1000;
        prefix++;
    }
    print_key("memory size", "%g %cB", size_bytes, *prefix);
}
