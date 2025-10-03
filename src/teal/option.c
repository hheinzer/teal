#include "option.h"

#include <getopt.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "assert.h"
#include "utils.h"

Option option = {0};

void option_init(int *argc, char ***argv)
{
    opterr = 0;
    int opt;
    while ((opt = getopt(*argc, *argv, "hqc:r:")) != -1) {
        switch (opt) {
            case '?':
            case 'h':
                print("usage: %s [-h] [-q] [-c capacity] [-r restart] ...\n", (*argv)[0]);
                MPI_Finalize();
                exit(EXIT_SUCCESS);
            case 'q': option.quiet = true; break;
            case 'c': option.capacity = str2size(optarg); break;
            case 'r': strcpy(option.restart.buf, optarg); break;
            default: assert(false);
        }
    }

    long num = 1;
    for (long i = optind; i < *argc; i++) {
        (*argv)[num++] = (*argv)[i];
    }
    *argc = num;
}
