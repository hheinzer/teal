#include <assert.h>

#include "parse.h"
#include "sync.h"
#include "teal.h"

Parse *parse_init(const char *fname)
{
    assert(fname);
    Parse *file = teal_alloc(1, sizeof(*file));
    int ret = MPI_File_open(sync.comm, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &file->handle);
    if (ret != MPI_SUCCESS) {
        teal_error("could not open file (%s)", fname);
    }
    MPI_File_set_errhandler(file->handle, MPI_ERRORS_ARE_FATAL);
    return file;
}
