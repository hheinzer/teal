#include <time.h>

#include "arena.h"
#include "sync.h"
#include "teal.h"
#include "utils.h"

void teal_finalize(void)
{
    char now[128];
    strftime(now, sizeof(now), "%a %b %e %T %Y", localtime(&(time_t){time(0)}));

    char min[128];
    size_to_str(min, sync_lmin(arena_size()));

    char max[128];
    size_to_str(max, sync_lmax(arena_size_max()));

    char sum[128];
    size_to_str(sum, sync_lsum(arena_size()));

    println("Goodbye, World!");
    println("\t stop time:              %s", now);
    println("\t arena min/max/sum size: %s / %s / %s", min, max, sum);

    arena_finalize();
    sync_finalize();
}
