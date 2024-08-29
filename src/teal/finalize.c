#include <stdio.h>
#include <time.h>

#include "option.h"
#include "print.h"
#include "sync.h"
#include "teal.h"

void teal_finalize(void)
{
    String s;
    strftime(s, sizeof(s), "%a %b %e %T %Y", localtime((time_t[]){time(0)}));

    if (sync.rank == 0 && !option.quiet) {
        printf("Goodbye, World!\n");
        print_key("end time", "%s", s);
    }

    sync_finalize();
}
