#include <time.h>

#include "arena.h"
#include "sync.h"
#include "teal.h"
#include "utils.h"

void teal_finalize(void)
{
    strbuf min = size2str(sync_lmin(arena_size()));
    strbuf max = size2str(sync_lmax(arena_size_max()));
    strbuf tot = size2str(sync_lsum(arena_size()));

    strbuf now;
    strftime(now.buf, sizeof(now), "%a %b %e %T %Y", localtime(&(time_t){time(0)}));

    print("Goodbye, World!\n");
    print("\t arena min/max/tot size: %s / %s / %s\n", min.buf, max.buf, tot.buf);
    print("\t stop time:              %s\n", now.buf);

    arena_finalize();
    sync_finalize();
}
