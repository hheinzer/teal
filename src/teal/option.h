#pragma once

#include <stdbool.h>

#include "teal.h"

typedef struct {
    bool quiet;
    bool verbose;
    int capacity;
    int num_refines;
    const char *restart;
} Option;

extern Option option;

void option_init(int *argc, char ***argv);
