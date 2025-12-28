#pragma once

extern struct Teal {
    int quiet;
} teal;

// Call exactly once at program begin.
void teal_init(int *argc, char ***argv);

// Call exactly once at program end.
void teal_deinit(void);

// Exit teal with a status code.
void teal_exit(int status) __attribute((noreturn));

// Print a formatted message from rank zero.
void teal_print(const char *fmt, ...) __attribute((format(printf, 1, 2)));

// Print a formatted error message and terminate.
void teal_error(const char *fmt, ...) __attribute((format(printf, 1, 2), noreturn));

// Allocate zeroed memory.
void *teal_alloc(long num, long size) __attribute((malloc));

// Reallocate or free an allocation.
void *teal_realloc(void *ptr, long num, long size);

// Free an allocation.
void teal_free(void *ptr);
