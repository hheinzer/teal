#include <time.h>

#include "arena.h"
#include "sync.h"
#include "teal.h"
#include "utils.h"

void teal_finalize(void)
{
    string min_size;
    size2str(min_size, sync_lmin(arena_size()));

    string max_size;
    size2str(max_size, sync_lmax(arena_size_max()));

    string tot_size;
    size2str(tot_size, sync_lsum(arena_size()));

    string now;
    strftime(now, sizeof(now), "%a %b %e %T %Y", localtime(&(time_t){time(0)}));

    print("Goodbye, World!\n");
    print("\t arena min/max/tot size: %s / %s / %s\n", min_size, max_size, tot_size);
    print("\t stop time:              %s\n", now);

    arena_finalize();
    sync_finalize();
}
