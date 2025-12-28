#include <assert.h>

#include "parse.h"
#include "teal.h"

void parse_deinit(Parse *file)
{
    assert(file);
    MPI_File_close(&file->handle);
    teal_free(file);
}
