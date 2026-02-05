#pragma once

#include <stddef.h>

extern struct Teal2 {
    int quiet;
    int verbose;
} teal2;

// Initialize teal; call once at program begin.
void teal2_init(int *argc, char ***argv);

// Finalize teal; call once at program end.
void teal2_deinit(void);

// Early exit teal with a status code.
void teal2_exit(int status) __attribute__((noreturn));

// Print a formatted message from a single rank to `stdout`.
void teal2_print(const char *fmt, ...) __attribute__((format(printf, 1, 2)));

// Print a formatted message from any rank to `stdout`.
void teal2_verbose(const char *fmt, ...) __attribute__((format(printf, 1, 2)));

// Print a formatted message from any rank to `stderr` and exit teal.
void teal2_error(const char *fmt, ...) __attribute__((format(printf, 1, 2), noreturn));

// Allocate memory with 64-byte alignment.
void *teal2_malloc(size_t size);

// Allocate `num * size` bytes and zero-initialize.
void *teal2_calloc(int num, int size) __attribute__((malloc));

// Free memory allocated with `teal2_malloc` or `teal2_calloc`.
void teal2_free(void *ptr);
