#pragma once

#include <stdbool.h>

#include "teal.h"

typedef struct {
    bool quiet;
    bool verbose;
    number capacity;
    number num_refines;
    const char *restart;
} Option;

extern Option option;

void option_init(int *argc, char ***argv);
