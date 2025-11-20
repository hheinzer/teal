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
    char max[128];
    char sum[128];

    println("Goodbye, World!");
    println("\t stop time              : %s", now);

    size_to_str(min, sync_lmin(arena_size()));
    size_to_str(max, sync_lmax(arena_size()));
    size_to_str(sum, sync_lsum(arena_size()));
    println("\t arena min/max/sum size : %s / %s / %s", min, max, sum);

    size_to_str(min, sync_lmin(arena_peak()));
    size_to_str(max, sync_lmax(arena_peak()));
    size_to_str(sum, sync_lsum(arena_peak()));
    println("\t arena min/max/sum peak : %s / %s / %s", min, max, sum);

    arena_finalize();
    sync_finalize();
}
