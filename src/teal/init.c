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

    string now;
    strftime(now, sizeof(now), "%a %b %e %T %Y", localtime(&(time_t){time(0)}));

    long capacity = option.capacity ? option.capacity : str2size("1G");
    arena_init(capacity);

    string cap;
    size2str(cap, capacity);

    print("Hello, World! This is teal!\n");
    print("\t start time:      %s\n", now);
    print("\t number of ranks: %d\n", sync.size);
    print("\t arena capacity:  %s\n", cap);
    print("\t restart file:    %s\n", *option.restart ? option.restart : "-");
}
