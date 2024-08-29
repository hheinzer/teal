#include "option.h"

#include <getopt.h>
#include <stdlib.h>

struct Option option = {0};

static void remove_read_arguments(int *argc, char ***argv);

void option_init(int *argc, char ***argv)
{
    int opt = 0;
    while ((opt = getopt(*argc, *argv, "qvr:")) != -1) {
        switch (opt) {
            case 'q': option.quiet = 1; break;
            case 'v': option.verbose = 1; break;
            case 'r': option.restart = optarg; break;
            default: abort();
        }
    }
    remove_read_arguments(argc, argv);
}

static void remove_read_arguments(int *argc, char ***argv)
{
    int n = 1;
    for (int i = optind; i < *argc; ++i) (*argv)[n++] = (*argv)[i];
    *argc = n;
}
