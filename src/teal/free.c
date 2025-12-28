#include <stdlib.h>

#include "teal.h"

void teal_free(void *ptr)
{
    free(ptr);
}
