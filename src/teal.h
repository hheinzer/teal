#pragma once

#include <stddef.h>

#define M_PI 3.14159265358979323846 /* pi */

extern struct teal {
    int quiet;
    int verbose;
    int partitioned;
} teal;

// Initialize teal; call once at program begin.
void teal_init(int *argc, char ***argv);

// Finalize teal; call once at program end.
void teal_deinit(void);

// Early exit teal with a status code.
void teal_exit(int status) __attribute__((noreturn));

// Print a formatted message from a single rank to `stdout`.
void teal_print(const char *fmt, ...) __attribute__((format(printf, 1, 2)));

// Print a formatted message from any rank to `stdout`.
void teal_verbose(const char *fmt, ...) __attribute__((format(printf, 1, 2)));

// Print a formatted message from any rank to `stderr` and exit teal.
void teal_error(const char *fmt, ...) __attribute__((format(printf, 1, 2), noreturn));

// Allocate `num * size` bytes and zero-initialize.
void *teal_calloc(int num, size_t size) __attribute__((malloc));

// Reallocate `num * size` bytes, preserving existing contents.
void *teal_realloc(void *ptr, int num, size_t size);

// Free memory allocated with `teal_calloc` or `teal_realloc`.
void teal_free(void *ptr);
