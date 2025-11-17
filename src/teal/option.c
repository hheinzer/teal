#include "option.h"

#include <getopt.h>
#include <mpi.h>
#include <stdlib.h>

#include "assert.h"
#include "sync.h"
#include "utils.h"

Option option = {0};

static void help(char ***argv, int status)
{
    println(
        "usage: %s"
        " [-h]"
        " [-q]"
        " [-v]"
        " [(-c | -C) capacity]"
        " [-n num_refines]"
        " [-r restart]"
        " ...",
        (*argv)[0]);
    sync_exit(status);
}

void option_init(int *argc, char ***argv)
{
    opterr = 0;
    int opt;
    while ((opt = getopt(*argc, *argv, ":hqvc:C:n:r:")) != -1) {
        switch (opt) {
            case 'h': help(argv, EXIT_SUCCESS); break;
            case 'q': option.quiet = true; break;
            case 'v': option.verbose = true; break;
            case 'c': option.capacity = str_to_size(optarg); break;
            case 'C': option.capacity = str_to_size(optarg) / sync.size; break;
            case 'n': option.num_refines = strtol(optarg, 0, 10); break;
            case 'r': option.restart = optarg; break;
            case '?':
                println("invalid option -- '%c'", optopt);
                help(argv, EXIT_FAILURE);
                break;
            case ':':
                println("option requires argument -- '%c'", optopt);
                help(argv, EXIT_FAILURE);
                break;
            default: assert(false);
        }
    }

    int num = 1;
    for (int i = optind; i < *argc; i++) {
        (*argv)[num++] = (*argv)[i];
    }
    *argc = num;
    (*argv)[num] = 0;
}
