#pragma once

#include <stdbool.h>

typedef struct {
    bool quiet;
    bool verbose;
    long capacity;
    long num_refines;
    char *restart;
} Option;

extern Option option;

void option_init(int *argc, char ***argv);
