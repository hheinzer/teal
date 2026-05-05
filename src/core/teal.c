#include "teal.h"

#include <assert.h>
#include <getopt.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "sanitizer.h"
#include "sync.h"
#include "utils.h"

_Static_assert(sizeof(int) == 4, "teal requires 32-bit int");
_Static_assert(sizeof(long) == 8, "teal requires 64-bit long");
_Static_assert(sizeof(double) == 8, "teal requires 64-bit double");

enum { ALIGN = 64 };

struct teal teal = {0};

static void print_help(char **argv)
{
    teal_print(
        "usage: %s [options] ...\n\n"
        "options:\n"
        "  %-20s show this help message and exit\n"
        "  %-20s disable normal and verbose printing (errors are still shown)\n"
        "  %-20s enable verbose printing\n"
        "  %-20s enable partitioned output files\n"
        "  %-20s restart simulation from output file\n",
        argv[0], "-h", "-q", "-v", "-p", "-r <file>");
}

static void parse_options(int *argc, char ***argv)
{
    opterr = 0;

    int opt;
    while ((opt = getopt(*argc, *argv, "hqvpr:")) != -1) {
        switch (opt) {
            case 'h': print_help(*argv); teal_exit(EXIT_SUCCESS);
            case 'q': teal.quiet = 1; break;
            case 'v': teal.verbose = 1; break;
            case 'p': teal.partitioned = 1; break;
            case 'r': teal.restart = optarg; break;
            default: print_help(*argv); teal_exit(EXIT_FAILURE);
        }
    }

    int num = 1;
    for (int i = optind; i < *argc; i++) {
        (*argv)[num++] = (*argv)[i];
    }

    *argc = num;
    (*argv)[num] = 0;
}

void teal_init(int *argc, char ***argv)
{
    assert(argc && argv);

    sync_init(argc, argv);
    parse_options(argc, argv);

    String now;
    strftime(now, sizeof(now), "%a %b %e %T %Y", localtime(&(time_t){time(0)}));

    teal_print("Hello, World! This is teal!");
    teal_print("\t start time      : %s", now);
    teal_print("\t number of ranks : %d", sync.size);
}

void teal_deinit(void)
{
    String now;
    strftime(now, sizeof(now), "%a %b %e %T %Y", localtime(&(time_t){time(0)}));

    teal_print("Goodbye, World!");
    teal_print("\t stop time : %s", now);

    sync_deinit();
}

void teal_exit(int status)
{
    sync_deinit();
    exit(status);
}

static void println(FILE *stream, const char *prefix, const char *fmt, va_list args)
{
    char beg[4 << 10];
    char *end = beg;

    if (prefix) {
        end += sprintf(end, "[%d] %s: ", sync.rank, prefix);
    }

    end += vsprintf(end, fmt, args);

    *end++ = '\n';
    *end = 0;

    fputs(beg, stream);
    fflush(stream);
}

void teal_print(const char *fmt, ...)
{
    assert(fmt);
    if (sync.rank == 0 && !teal.quiet) {
        va_list args;
        va_start(args, fmt);
        println(stdout, 0, fmt, args);
        va_end(args);
    }
}

void teal_verbose(const char *fmt, ...)
{
    assert(fmt);
    if (!teal.quiet && teal.verbose) {
        va_list args;
        va_start(args, fmt);
        println(stdout, "VERBOSE", fmt, args);
        va_end(args);
    }
}

void teal_error(const char *fmt, ...)
{
    assert(fmt);

    va_list args;
    va_start(args, fmt);
    println(stderr, "ERROR", fmt, args);
    va_end(args);

    teal_exit(EXIT_FAILURE);
}

typedef struct {
    void *base;
    size_t bytes;
} Alloc;

static void *teal_malloc(int num, size_t size)
{
    size_t extra = sizeof(Alloc) + (ALIGN - 1);
    if ((size_t)num > (SIZE_MAX - extra) / size) {
        teal_error("overflow (%d, %zu)", num, size);
    }

    size_t bytes = num * size;

    char *base = malloc(bytes + extra);
    if (!base) {
        teal_error("malloc failure (%zu)", bytes + extra);
    }

    char *beg = base + sizeof(Alloc);
    char *ptr = beg + (-(uintptr_t)beg & (ALIGN - 1));

    ((Alloc *)ptr)[-1].base = base;
    ((Alloc *)ptr)[-1].bytes = bytes;

    MAKE_REGION_NOACCESS(base, ptr - base);

    return ptr;
}

void *teal_calloc(int num, size_t size)
{
    assert(num >= 0 && size > 0);

    if (num == 0) {
        return 0;
    }

    void *ptr = teal_malloc(num, size);

    return memset(ptr, 0, num * size);
}

void *teal_realloc(void *ptr, int num, size_t size)
{
    assert(num >= 0 && size > 0);

    if (num == 0) {
        teal_free(ptr);
        return 0;
    }

    if (!ptr) {
        return teal_malloc(num, size);
    }

    MAKE_REGION_ADDRESSABLE((char *)ptr - sizeof(Alloc), sizeof(Alloc));

    void *new = teal_malloc(num, size);

    size_t bytes = num * size;
    if (((Alloc *)ptr)[-1].bytes < bytes) {
        bytes = ((Alloc *)ptr)[-1].bytes;
    }

    memcpy(new, ptr, bytes);

    free(((Alloc *)ptr)[-1].base);

    return new;
}

void teal_free(void *ptr)
{
    if (!ptr) {
        return;
    }

    MAKE_REGION_ADDRESSABLE((char *)ptr - sizeof(Alloc), sizeof(Alloc));

    free(((Alloc *)ptr)[-1].base);
}
