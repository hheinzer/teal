#include <stdio.h>
#include <time.h>

#include "option.h"
#include "print.h"
#include "sync.h"
#include "teal.h"

void teal_initialize(int *argc, char ***argv)
{
    sync_init(argc, argv);
    option_init(argc, argv);

    String s;
    strftime(s, sizeof(s), "%a %b %e %T %Y", localtime((time_t[]){time(0)}));

    if (sync.rank == 0 && !option.quiet) {
        printf("Hello, World! This is teal!\n");
        print_key("start time", "%s", s);
        print_key("number of ranks", "%d", sync.size);
    }
}
