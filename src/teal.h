#pragma once

#include <stdbool.h>

typedef double scalar;

typedef struct {
    scalar x, y, z;
} vector;

typedef struct {
    vector x, y, z;
} matrix;

typedef struct {
    long x, y, z;
} tuple;

// Call exactly once at program start, before any teal API.
void teal_initialize(int *argc, char ***argv);

// Call exactly once at program end, after teal usage.
void teal_finalize(void);

// Safely exit teal with a status code.
void teal_exit(int status) __attribute((noreturn));
