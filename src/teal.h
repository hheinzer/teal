#pragma once

#include <stdbool.h>
#include <stddef.h>

typedef ptrdiff_t number;

typedef double scalar;

typedef struct {
    scalar x, y, z;
} vector;

typedef struct {
    vector x, y, z;
} matrix;

typedef struct {
    number x, y, z;
} tuple;

typedef struct {
    bool x, y, z;
} flags;

/* Call exactly once at program start, before any teal API. */
void teal_init(int *argc, char ***argv);

/* Call exactly once at program end, after teal usage. */
void teal_finalize(void);
