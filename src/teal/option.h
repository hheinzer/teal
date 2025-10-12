#pragma once

#include "teal.h"

typedef struct {
    bool quiet;
    long capacity;
    long num_refine;
    strbuf restart;
} Option;

extern Option option;

void option_init(int *argc, char ***argv);
