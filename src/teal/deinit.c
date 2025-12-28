#include <time.h>

#include "sync.h"
#include "teal.h"

void teal_deinit(void)
{
    char now[128];
    strftime(now, sizeof(now), "%a %b %e %T %Y", localtime(&(time_t){time(0)}));

    teal_print("Goodbye, World!");
    teal_print("\t stop time : %s", now);

    sync_deinit();
}
