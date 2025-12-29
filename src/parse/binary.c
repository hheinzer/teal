#include <assert.h>
#include <mpi.h>

#include "parse.h"
#include "sync.h"
#include "teal.h"

int parse_binary(Parse *file, void *buf, int num, MPI_Datatype datatype)
{
    assert(file && (buf || num == 0) && num >= 0);
    if (num == 0) {
        return 0;
    }
    int count = 0;
    if (sync.rank == 0) {
        MPI_Status status;
        MPI_File_read_at(file->handle, file->offset, buf, num, datatype, &status);
        MPI_Get_count(&status, datatype, &count);
        if (count <= 0) {
            teal_error("invalid read (probably end-of-file)");
        }
        int size = 0;
        MPI_Type_size(datatype, &size);
        if (size <= 0) {
            teal_error("invalid type size (%d)", size);
        }
        file->offset += (MPI_Offset)count * size;
    }
    MPI_Bcast(&count, 1, MPI_INT, 0, sync.comm);
    MPI_Bcast(buf, count, datatype, 0, sync.comm);
    MPI_Bcast(&file->offset, 1, MPI_OFFSET, 0, sync.comm);
    return count;
}

int parse_binary_split(Parse *file, void *buf, int num, MPI_Datatype datatype)
{
    assert(file && (buf || num == 0) && num >= 0);
    int size = 0;
    MPI_Type_size(datatype, &size);
    if (size <= 0) {
        teal_error("invalid type size (%d)", size);
    }

    long prefix = sync_exsum(num);
    MPI_Offset offset = file->offset + (prefix * size);

    int count = 0;
    if (num > 0) {
        MPI_Status status;
        MPI_File_read_at(file->handle, offset, buf, num, datatype, &status);
        MPI_Get_count(&status, datatype, &count);
        if (count <= 0) {
            teal_error("invalid read (probably end-of-file)");
        }
    }

    long total = sync_all_lsum(count);
    file->offset += total * size;
    return count;
}
