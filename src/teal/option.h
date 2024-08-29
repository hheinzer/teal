#pragma once

extern struct Option {
    int quiet, verbose;
    const char *restart;
} option;

void option_init(int *argc, char ***argv);
