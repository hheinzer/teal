#include <time.h>

#include "arena.h"
#include "option.h"
#include "sync.h"
#include "teal.h"
#include "utils.h"

void teal_init(int *argc, char ***argv)
{
    sync_init(argc, argv);
    option_init(argc, argv);

    char now[128];
    strftime(now, sizeof(now), "%a %b %e %T %Y", localtime(&(time_t){time(0)}));

    long capacity = option.capacity ? option.capacity : str_to_size("1G");
    arena_init(capacity);

    char cap[128];
    size_to_str(cap, capacity);

    println("Hello, World! This is teal!");
    println("\t start time:        %s", now);
    println("\t number of ranks:   %d", sync.size);
    println("\t arena capacity:    %s", cap);
    println("\t number of refines: %ld", option.num_refines);
    println("\t restart file:      %s", option.restart);
}
