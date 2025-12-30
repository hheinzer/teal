#include <assert.h>
#include <byteswap.h>
#include <limits.h>
#include <mpi.h>
#include <stdint.h>

#include "parse.h"
#include "sync.h"
#include "teal.h"

// Swap byte order for fixed-size elements.
static void swap_bytes(void *buf, int num, int size)
{
    switch (size) {
        case 1: break;
        case 2: {
            uint16_t *u16 = buf;
            for (long i = 0; i < num; i++) {
                u16[i] = bswap_16(u16[i]);
            }
            break;
        }
        case 4: {
            uint32_t *u32 = buf;
            for (long i = 0; i < num; i++) {
                u32[i] = bswap_32(u32[i]);
            }
            break;
        }
        case 8: {
            uint64_t *u64 = buf;
            for (long i = 0; i < num; i++) {
                u64[i] = bswap_64(u64[i]);
            }
            break;
        }
        default: teal_error("invalid size (%d)", size);
    }
}

int parse_binary(Parse *file, void *buf, int num, MPI_Datatype datatype, int mode)
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
        if (mode & SWAP) {
            swap_bytes(buf, count, size);
        }
        file->offset += (MPI_Offset)count * size;
    }
    MPI_Bcast(&count, 1, MPI_INT, 0, sync.comm);
    MPI_Bcast(buf, count, datatype, 0, sync.comm);
    MPI_Bcast(&file->offset, 1, MPI_OFFSET, 0, sync.comm);
    return count;
}

int parse_binary_split(Parse *file, void *buf, int num, MPI_Datatype datatype, int mode)
{
    assert(file && (buf || num == 0) && num >= 0);

    int size = 0;
    MPI_Type_size(datatype, &size);
    if (size <= 0) {
        teal_error("invalid type size (%d)", size);
    }

    long prefix = sync_lexsum(num);
    MPI_Offset offset = file->offset + (prefix * size);

    int count = 0;
    if (num > 0) {
        MPI_Status status;
        MPI_File_read_at(file->handle, offset, buf, num, datatype, &status);
        MPI_Get_count(&status, datatype, &count);
        if (count <= 0) {
            teal_error("invalid read (probably end-of-file)");
        }
        if (mode & SWAP) {
            swap_bytes(buf, count, size);
        }
    }

    long total = sync_lsum(count);
    file->offset += total * size;
    return count;
}
