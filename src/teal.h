#pragma once

#include <stddef.h>

extern struct Teal {
    int quiet;
} teal;

// Initialize teal and MPI once at program start.
void teal_init(int *argc, char ***argv);

// Finalize teal and MPI once at program end.
void teal_deinit(void);

// Exit teal with a status code.
void teal_exit(int status) __attribute((noreturn));

// Print a formatted message from rank zero.
void teal_print(const char *fmt, ...) __attribute((format(printf, 1, 2)));

// Print a formatted error message and terminate.
void teal_error(const char *fmt, ...) __attribute((format(printf, 1, 2), noreturn));

// Allocate zeroed memory for num elements.
void *teal_alloc(int num, size_t size) __attribute((malloc));

// Resize an allocation or free it when num is 0.
void *teal_realloc(void *ptr, int num, size_t size);

// Free an allocation.
void teal_free(void *ptr);
