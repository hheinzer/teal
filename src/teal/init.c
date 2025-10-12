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

    strbuf now;
    strftime(now.buf, sizeof(now), "%a %b %e %T %Y", localtime(&(time_t){time(0)}));

    long capacity = option.capacity ? option.capacity : str2size("1G");
    arena_init(capacity);

    print("Hello, World! This is teal!\n");
    print("\t start time:        %s\n", now.buf);
    print("\t number of ranks:   %d\n", sync.size);
    print("\t arena capacity:    %s\n", size2str(capacity).buf);
    print("\t number of refines: %ld\n", option.num_refines);
    print("\t restart file:      %s\n", *option.restart.buf ? option.restart.buf : "-");
}
