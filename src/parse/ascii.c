#include <assert.h>
#include <errno.h>
#include <mpi.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "parse.h"
#include "space.h"
#include "sync.h"
#include "teal.h"

enum { SIZE = 4096 };

typedef struct {
    char data[SIZE + 1];
    char *beg;
    char *end;
    MPI_Offset offset;
    int eof;
} Buffer;

// Refill the buffer and return bytes read or 0 on EOF.
static int refill(Buffer *buffer, Parse *file)
{
    if (buffer->eof) {
        return 0;
    }
    if (buffer->beg > buffer->data) {
        // slide unread data to the front of the buffer
        long consumed = buffer->beg - buffer->data;
        long remaining = buffer->end - buffer->beg;
        memmove(buffer->data, buffer->beg, (size_t)remaining);
        buffer->offset += consumed;
        buffer->beg = buffer->data;
        buffer->end = buffer->data + remaining;
    }
    int used = (int)(buffer->end - buffer->data);
    int available = SIZE - used;
    if (available <= 0) {
        return 0;
    }
    MPI_Status status;
    MPI_File_read_at(file->handle, buffer->offset + used, buffer->end, available, MPI_CHAR,
                     &status);
    int count = 0;
    MPI_Get_count(&status, MPI_CHAR, &count);
    if (count <= 0) {
        buffer->eof = 1;
        return 0;
    }
    buffer->end += count;
    return count;
}

// Return the next token in place or 0 on EOF.
static char *next(Buffer *buffer, Parse *file)
{
    if (buffer->beg >= buffer->end) {
        if (!refill(buffer, file)) {
            return 0;
        }
    }
    while (1) {
        buffer->beg = skip_space(buffer->beg, buffer->end);
        if (buffer->beg < buffer->end) {
            break;
        }
        if (!refill(buffer, file)) {
            return 0;
        }
    }
    char *tok = buffer->beg;
    while (1) {
        // scan token until a delimiter or eof is found
        char *end = find_space(tok, buffer->end);
        if (end < buffer->end) {
            *end = 0;
            buffer->beg = end + 1;
            return tok;
        }

        if (buffer->eof) {
            *buffer->end = 0;
            buffer->beg = buffer->end;
            return tok;
        }

        // refill while keeping the current token start valid
        buffer->beg = tok;
        if (buffer->beg == buffer->data && buffer->end >= buffer->data + SIZE) {
            teal_error("token too long (exceeds %d bytes)", SIZE);
        }
        if (!refill(buffer, file)) {
            if (!buffer->eof) {
                teal_error("token too long (exceeds %d bytes)", SIZE);
            }
            continue;
        }
        tok = buffer->beg;
    }
}

// Convert one token into the requested MPI datatype.
static void convert(const char *token, void *buf, int index, MPI_Datatype datatype)
{
    errno = 0;
    char *end = 0;
    if (datatype == MPI_INT8_T) {
        ((int8_t *)buf)[index] = (int8_t)strtol(token, &end, 10);
    }
    else if (datatype == MPI_INT16_T) {
        ((int16_t *)buf)[index] = (int16_t)strtol(token, &end, 10);
    }
    else if (datatype == MPI_INT32_T) {
        ((int32_t *)buf)[index] = (int32_t)strtol(token, &end, 10);
    }
    else if (datatype == MPI_INT64_T) {
        ((int64_t *)buf)[index] = (int64_t)strtoll(token, &end, 10);
    }
    else if (datatype == MPI_UINT8_T) {
        ((uint8_t *)buf)[index] = (uint8_t)strtoul(token, &end, 10);
    }
    else if (datatype == MPI_UINT16_T) {
        ((uint16_t *)buf)[index] = (uint16_t)strtoul(token, &end, 10);
    }
    else if (datatype == MPI_UINT32_T) {
        ((uint32_t *)buf)[index] = (uint32_t)strtoul(token, &end, 10);
    }
    else if (datatype == MPI_UINT64_T) {
        ((uint64_t *)buf)[index] = (uint64_t)strtoull(token, &end, 10);
    }
    else if (datatype == MPI_FLOAT) {
        ((float *)buf)[index] = strtof(token, &end);
    }
    else if (datatype == MPI_DOUBLE) {
        ((double *)buf)[index] = strtod(token, &end);
    }
    else {
        char name[MPI_MAX_OBJECT_NAME];
        int len = 0;
        MPI_Type_get_name(datatype, name, &len);
        teal_error("unsupported datatype (%.*s)", len, name);
    }
    if (errno != 0 || !end || end == token || *end != 0) {
        char name[MPI_MAX_OBJECT_NAME];
        int len = 0;
        MPI_Type_get_name(datatype, name, &len);
        teal_error("could not parse token (%s) as %.*s", token, len, name);
    }
}

int parse_ascii(Parse *file, void *buf, int num, MPI_Datatype datatype)
{
    assert(file && buf && num > 0);
    int count = 0;
    if (sync.rank == 0) {
        Buffer buffer = {
            .beg = buffer.data,
            .end = buffer.data,
            .offset = file->offset,
        };
        for (int i = 0; i < num; i++) {
            const char *tok = next(&buffer, file);
            if (!tok) {
                break;
            }
            convert(tok, buf, i, datatype);
            count += 1;
        }
        file->offset = buffer.offset + buffer.beg - buffer.data;
    }
    MPI_Bcast(&count, 1, MPI_INT, 0, sync.comm);
    MPI_Bcast(buf, count, datatype, 0, sync.comm);
    MPI_Bcast(&file->offset, 1, MPI_OFFSET, 0, sync.comm);
    return count;
}
