#include "option.h"

#include <getopt.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "assert.h"
#include "teal/sync.h"
#include "utils.h"

Option option = {0};

void option_init(int *argc, char ***argv)
{
    opterr = 0;
    int opt;
    while ((opt = getopt(*argc, *argv, "hqc:C:n:r:")) != -1) {
        switch (opt) {
            case '?':
            case 'h':
                print(
                    "usage: %s"
                    " [-h]"
                    " [-q]"
                    " [-c | -C capacity]"
                    " [-n num_refine]"
                    " [-r restart]"
                    " ...\n",
                    (*argv)[0]);
                MPI_Finalize();
                exit(EXIT_SUCCESS);
            case 'q': option.quiet = true; break;
            case 'c': option.capacity = str2size(optarg); break;
            case 'C': option.capacity = str2size(optarg) / sync.size; break;
            case 'n': option.num_refine = strtol(optarg, 0, 0); break;
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
