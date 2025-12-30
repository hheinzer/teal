#include <assert.h>
#include <getopt.h>
#include <stdlib.h>
#include <time.h>

#include "sync.h"
#include "teal.h"

struct Teal teal = {0};

// Print usage and terminate with the given status.
static void print_help(char **argv, int status)
{
    teal_print(
        "usage: %s"
        " [-h]"
        " [-q]"
        " ...",
        argv[0]);
    teal_exit(status);
}

// Parse teal command line options.
static void parse_options(int argc, char **argv)
{
    opterr = 0;

    int opt;
    while ((opt = getopt(argc, argv, "hq")) != -1) {
        switch (opt) {
            case 'h': print_help(argv, EXIT_SUCCESS); break;
            case 'q': teal.quiet = 1; break;
            default: print_help(argv, EXIT_FAILURE);
        }
    }
}

// Remove parsed options from argv and update argc.
static void remove_parsed_options(int *argc, char ***argv)
{
    int num = 1;
    for (int i = optind; i < *argc; i++) {
        (*argv)[num++] = (*argv)[i];
    }
    *argc = num;
    (*argv)[num] = 0;
}

void teal_init(int *argc, char ***argv)
{
    assert(argc && argv);

    sync_init(argc, argv);

    if (sizeof(int) != 4) {
        teal_error("unexpected sizeof(int) (%zu)", sizeof(int));
    }
    if (sizeof(long) != 8) {
        teal_error("unexpected sizeof(long) (%zu)", sizeof(long));
    }

    parse_options(*argc, *argv);
    remove_parsed_options(argc, argv);

    char now[128];
    strftime(now, sizeof(now), "%a %b %e %T %Y", localtime(&(time_t){time(0)}));

    teal_print("Hello, World! This is teal!");
    teal_print("\t start time      : %s", now);
    teal_print("\t number of ranks : %d", sync.size);
}
