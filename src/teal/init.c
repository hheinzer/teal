#include <time.h>

#include "arena.h"
#include "option.h"
#include "sync.h"
#include "teal.h"
#include "utils.h"

#ifndef COMPILER
#define COMPILER "unknown"
#endif

#ifndef COMMIT
#define COMMIT "unknown"
#endif

#ifndef CONFIG
#define CONFIG "unknown"
#endif

void teal_init(int *argc, char ***argv)
{
    sync_init(argc, argv);
    option_init(argc, argv);

    char now[128];
    strftime(now, sizeof(now), "%a %b %e %T %Y", localtime(&(time_t){time(0)}));

    number capacity = option.capacity ? option.capacity : str_to_size("1G");
    arena_init(capacity);

    char cap[128];
    size_to_str(cap, capacity);

    println("Hello, World! This is teal!");
    println("\t start time        : %s", now);
    println("\t compiler          : %s", COMPILER);
    println("\t commit            : %s", COMMIT);
    println("\t config            : %s", CONFIG);
    println("\t number of ranks   : %d", sync.size);
    println("\t arena capacity    : %s", cap);
    if (option.num_refines > 0) {
        println("\t number of refines : %td", option.num_refines);
    }
    if (option.restart) {
        println("\t restart file      : %s", option.restart);
    }
}
