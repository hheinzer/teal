#include <assert.h>
#include <getopt.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "sanitizer.h"
#include "sync2.h"
#include "teal2.h"

struct Teal2 teal2 = {0};

static void print_help(char **argv)
{
    teal2_print("usage: %s [-h] [-q] [-v] ...", argv[0]);
}

static void parse_options(int *argc, char ***argv)
{
    opterr = 0;

    int opt;
    while ((opt = getopt(*argc, *argv, "hqv")) != -1) {
        switch (opt) {
            case 'h': print_help(*argv); teal2_exit(EXIT_SUCCESS);
            case 'q': teal2.quiet = 1; break;
            case 'v': teal2.verbose = 1; break;
            default: print_help(*argv); teal2_exit(EXIT_FAILURE);
        }
    }

    int num = 1;
    for (int i = optind; i < *argc; i++) {
        (*argv)[num++] = (*argv)[i];
    }

    *argc = num;
    (*argv)[num] = 0;
}

void teal2_init(int *argc, char ***argv)
{
    assert(argc && argv);

    sync2_init(argc, argv);
    parse_options(argc, argv);

    char now[128];
    strftime(now, sizeof(now), "%a %b %e %T %Y", localtime(&(time_t){time(0)}));

    teal2_print("Hello, World! This is teal!");
    teal2_print("\t start time      : %s", now);
    teal2_print("\t number of ranks : %d", sync2.size);
}

void teal2_deinit(void)
{
    char now[128];
    strftime(now, sizeof(now), "%a %b %e %T %Y", localtime(&(time_t){time(0)}));

    teal2_print("Goodbye, World!");
    teal2_print("\t stop time : %s", now);

    sync2_deinit();
}

void teal2_exit(int status)
{
    sync2_deinit();
    exit(status);
}

static void println(FILE *stream, const char *prefix, const char *fmt, va_list args)
{
    char beg[4 << 10];
    char *end = beg;
    if (prefix) {
        end += sprintf(end, "[%d] %s: ", sync2.rank, prefix);
    }
    end += vsprintf(end, fmt, args);
    *end++ = '\n';
    *end = 0;
    fputs(beg, stream);
    fflush(stream);
}

void teal2_print(const char *fmt, ...)
{
    assert(fmt);
    if (sync2.rank == 0 && !teal2.quiet) {
        va_list args;
        va_start(args, fmt);
        println(stdout, 0, fmt, args);
        va_end(args);
    }
}

void teal2_verbose(const char *fmt, ...)
{
    assert(fmt);
    if (!teal2.quiet && teal2.verbose) {
        va_list args;
        va_start(args, fmt);
        println(stdout, "VERBOSE", fmt, args);
        va_end(args);
    }
}

void teal2_error(const char *fmt, ...)
{
    assert(fmt);
    va_list args;
    va_start(args, fmt);
    println(stdout, "ERROR", fmt, args);
    va_end(args);
    teal2_exit(EXIT_FAILURE);
}

typedef struct {
    void *base;
    size_t bytes;
} Alloc;

void *teal2_malloc(size_t bytes)
{
    if (bytes == 0) {
        return 0;
    }

    enum { ALIGN = 64 };

    size_t extra = sizeof(Alloc) + (ALIGN - 1);
    if (bytes > SIZE_MAX - extra) {
        teal2_error("overflow (%zu)", bytes);
    }

    char *base = malloc(bytes + extra);
    if (!base) {
        teal2_error("malloc failure (%zu)", bytes);
    }

    char *beg = base + sizeof(Alloc);
    char *ptr = beg + (-(uintptr_t)beg & (ALIGN - 1));

    ((Alloc *)ptr)[-1].base = base;
    ((Alloc *)ptr)[-1].bytes = bytes;

    MAKE_REGION_NOACCESS(base, ptr - base);

    return ptr;
}

void *teal2_calloc(int num, int size)
{
    assert(num >= 0 && size > 0);

    if (num == 0) {
        return 0;
    }

    if ((size_t)num > SIZE_MAX / size) {
        teal2_error("overflow (%d, %d)", num, size);
    }

    size_t bytes = (size_t)num * size;
    void *ptr = teal2_malloc(bytes);

    return memset(ptr, 0, bytes);
}

void teal2_free(void *ptr)
{
    if (!ptr) {
        return;
    }

    MAKE_REGION_ADDRESSABLE((char *)ptr - sizeof(Alloc), sizeof(Alloc));

    free(((Alloc *)ptr)[-1].base);
}

void *teal2_recalloc(void *ptr, int num, int size)
{
    assert(num >= 0 && size > 0);

    if (!ptr) {
        return teal2_calloc(num, size);
    }

    if (num == 0) {
        teal2_free(ptr);
        return 0;
    }

    void *new = teal2_calloc(num, size);
    size_t bytes = (size_t)num * size;

    MAKE_REGION_ADDRESSABLE((char *)ptr - sizeof(Alloc), sizeof(Alloc));

    if (((Alloc *)ptr)[-1].bytes < bytes) {
        bytes = ((Alloc *)ptr)[-1].bytes;
    }

    memcpy(new, ptr, bytes);

    free(((Alloc *)ptr)[-1].base);

    return new;
}
