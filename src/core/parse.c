#include "parse.h"

#include <assert.h>
#include <byteswap.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "sync.h"
#include "teal.h"

enum { SIZE = 4 << 10 };

MPI_File parse_open(const char *fname)
{
    assert(fname);

    MPI_File file;

    if (MPI_File_open(sync.comm, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &file)) {
        teal_error("could not open file (%s)", fname);
    }

    if (MPI_File_set_errhandler(file, MPI_ERRORS_ARE_FATAL)) {
        teal_error("could not install error handler");
    }

    return file;
}

void parse_close(MPI_File file)
{
    MPI_File_close(&file);
}

static char *skip_space(char *beg, const char *end)
{
    while (beg < end && isspace((unsigned char)*beg)) {
        beg += 1;
    }
    return beg;
}

static int isquote(char chr)
{
    switch (chr) {
        case '"':
        case '\'':
        case '`': return 1;
        default: return 0;
    }
}

static char *find_quote(char *beg, const char *end)
{
    char quote = *beg;
    beg += 1;
    while (beg < end && *beg != quote) {
        beg += 1;
    }
    return beg;
}

static char *find_space(char *beg, const char *end)
{
    while (beg < end && !isspace((unsigned char)*beg)) {
        beg += 1;
    }
    return beg;
}

static int read_token(MPI_Offset *offset, char *str, int size, int count)
{
    char *end = str + count;
    char *beg = skip_space(str, end);
    if (beg >= end) {
        teal_error("invalid string (only spaces)");
    }

    if (isquote(*beg)) {
        char *quote = find_quote(beg, end);
        if (quote >= end) {
            teal_error("unterminated quoted string");
        }
        end = quote;
        beg += 1;
    }
    else {
        char *space = find_space(beg, end);
        if (space >= end && count == size) {
            teal_error("missing delimiter (buffer too small)");
        }
        end = space;
    }

    *offset += end - str;
    if (end < str + count) {
        *offset += 1;
    }

    ptrdiff_t len = end - beg;
    memmove(str, beg, len);
    memset(str + len, 0, size - len);

    assert(len <= INT_MAX);
    return (int)len;
}

int parse_string(MPI_File file, char *str, int size)
{
    assert(str && size > 0);

    int len = 0;
    MPI_Offset offset;
    if (sync.rank == 0) {
        MPI_File_get_position(file, &offset);

        MPI_Status status;
        MPI_File_read_at(file, offset, str, size, MPI_CHAR, &status);

        int count = 0;
        MPI_Get_count(&status, MPI_CHAR, &count);

        if (count == 0) {
            memset(str, 0, size);
        }
        else {
            len = read_token(&offset, str, size, count);
        }
    }

    MPI_Bcast(&len, 1, MPI_INT, 0, sync.comm);
    MPI_Bcast(str, size, MPI_CHAR, 0, sync.comm);

    MPI_Bcast(&offset, 1, MPI_OFFSET, 0, sync.comm);
    MPI_File_seek(file, offset, MPI_SEEK_SET);

    return len;
}

typedef struct {
    MPI_File file;
    MPI_Offset offset;
    char base[SIZE + 1];
    char *beg;
    char *end;
    int eof;
} Buffer;

static int refill(Buffer *buffer)
{
    if (buffer->eof) {
        return 0;
    }

    if (buffer->beg > buffer->base) {
        ptrdiff_t consumed = buffer->beg - buffer->base;
        ptrdiff_t remaining = buffer->end - buffer->beg;
        memmove(buffer->base, buffer->beg, remaining);
        buffer->offset += consumed;
        buffer->beg = buffer->base;
        buffer->end = buffer->base + remaining;
    }

    ptrdiff_t used = buffer->end - buffer->base;
    assert(used <= INT_MAX);

    int available = SIZE - (int)used;
    if (available <= 0) {
        return 0;
    }

    MPI_Status status;
    MPI_File_read_at(buffer->file, buffer->offset + used, buffer->end, available, MPI_CHAR,
                     &status);

    int count = 0;
    MPI_Get_count(&status, MPI_CHAR, &count);

    if (count == 0) {
        buffer->eof = 1;
        return 0;
    }

    buffer->end += count;

    return count;
}

static char *next_token(Buffer *buffer)
{
    if (buffer->beg >= buffer->end) {
        if (!refill(buffer)) {
            return 0;
        }
    }

    while (1) {
        buffer->beg = skip_space(buffer->beg, buffer->end);
        if (buffer->beg < buffer->end) {
            break;
        }
        if (!refill(buffer)) {
            return 0;
        }
    }

    char *token = buffer->beg;
    while (1) {
        char *end = find_space(token, buffer->end);
        if (end < buffer->end) {
            *end = 0;
            buffer->beg = end + 1;
            return token;
        }

        if (buffer->eof) {
            *buffer->end = 0;
            buffer->beg = buffer->end;
            return token;
        }

        buffer->beg = token;
        if (buffer->beg == buffer->base && buffer->end >= buffer->base + SIZE) {
            teal_error("token too long (exceeds %d bytes)", SIZE);
        }

        if (!refill(buffer)) {
            if (!buffer->eof) {
                teal_error("token too long (exceeds %d bytes)", SIZE);
            }
            continue;
        }

        token = buffer->beg;
    }
}

static long long str_to_signed(const char *token, char **end, long long min, long long max)
{
    long long val = strtoll(token, end, 10);
    if (val < min || max < val) {
        teal_error("value %s out of range [%lld, %lld]", token, min, max);
    }
    return val;
}

static unsigned long long str_to_unsigned(const char *token, char **end, unsigned long long max)
{
    unsigned long long val = strtoull(token, end, 10);
    if (max < val) {
        teal_error("value %s out of range [0, %llu]", token, max);
    }
    return val;
}

static void convert(const char *token, void *buf, int idx, MPI_Datatype type)
{
    errno = 0;
    char *end = 0;
    if (type == MPI_INT8_T) {
        ((int8_t *)buf)[idx] = (int8_t)str_to_signed(token, &end, INT8_MIN, INT8_MAX);
    }
    else if (type == MPI_INT16_T) {
        ((int16_t *)buf)[idx] = (int16_t)str_to_signed(token, &end, INT16_MIN, INT16_MAX);
    }
    else if (type == MPI_INT32_T) {
        ((int32_t *)buf)[idx] = (int32_t)str_to_signed(token, &end, INT32_MIN, INT32_MAX);
    }
    else if (type == MPI_INT64_T) {
        ((int64_t *)buf)[idx] = (int64_t)str_to_signed(token, &end, INT64_MIN, INT64_MAX);
    }
    else if (type == MPI_UINT8_T) {
        ((uint8_t *)buf)[idx] = (uint8_t)str_to_unsigned(token, &end, UINT8_MAX);
    }
    else if (type == MPI_UINT16_T) {
        ((uint16_t *)buf)[idx] = (uint16_t)str_to_unsigned(token, &end, UINT16_MAX);
    }
    else if (type == MPI_UINT32_T) {
        ((uint32_t *)buf)[idx] = (uint32_t)str_to_unsigned(token, &end, UINT32_MAX);
    }
    else if (type == MPI_UINT64_T) {
        ((uint64_t *)buf)[idx] = (uint64_t)str_to_unsigned(token, &end, UINT64_MAX);
    }
    else if (type == MPI_FLOAT) {
        ((float *)buf)[idx] = strtof(token, &end);
    }
    else if (type == MPI_DOUBLE) {
        ((double *)buf)[idx] = strtod(token, &end);
    }
    else {
        char name[MPI_MAX_OBJECT_NAME];
        int len = 0;
        MPI_Type_get_name(type, name, &len);
        teal_error("unsupported type (%.*s)", len, name);
    }

    if (errno != 0 || !end || end == token || *end != 0) {
        char name[MPI_MAX_OBJECT_NAME];
        int len = 0;
        MPI_Type_get_name(type, name, &len);
        teal_error("could not parse token (%s) as %.*s", token, len, name);
    }
}

static int read_tokens(Buffer *buffer, void *buf, int num, MPI_Datatype type)
{
    int count = 0;
    for (int i = 0; i < num; i++) {
        char *token = next_token(buffer);
        if (!token) {
            break;
        }
        convert(token, buf, i, type);
        count += 1;
    }
    return count;
}

void parse_ascii(MPI_File file, void *buf, int num, MPI_Datatype type)
{
    assert((buf || num == 0) && num >= 0);

    if (num == 0) {
        return;
    }

    MPI_Offset offset;
    if (sync.rank == 0) {
        MPI_File_get_position(file, &offset);

        Buffer buffer = {0};
        buffer.file = file;
        buffer.offset = offset;
        buffer.beg = buffer.base;
        buffer.end = buffer.base;

        int count = read_tokens(&buffer, buf, num, type);
        if (count != num) {
            teal_error("short read (%d / %d)", count, num);
        }

        offset = buffer.offset + buffer.beg - buffer.base;
    }

    MPI_Bcast(buf, num, type, 0, sync.comm);

    MPI_Bcast(&offset, 1, MPI_OFFSET, 0, sync.comm);
    MPI_File_seek(file, offset, MPI_SEEK_SET);
}

void parse_ascii_split(MPI_File file, void *buf, int num, MPI_Datatype type)
{
    assert((buf || num == 0) && num >= 0);

    for (int rank = 0; rank < sync.size; rank++) {
        MPI_Offset offset;
        if (sync.rank == rank) {
            MPI_File_get_position(file, &offset);

            if (num > 0) {
                Buffer buffer = {0};
                buffer.file = file;
                buffer.offset = offset;
                buffer.beg = buffer.base;
                buffer.end = buffer.base;

                int count = read_tokens(&buffer, buf, num, type);
                if (count != num) {
                    teal_error("short read (%d / %d)", count, num);
                }

                offset = buffer.offset + buffer.beg - buffer.base;
            }
        }
        MPI_Bcast(&offset, 1, MPI_OFFSET, rank, sync.comm);
        MPI_File_seek(file, offset, MPI_SEEK_SET);
    }
}

static void swap_bytes(void *buf, int num, int size)
{
    switch (size) {
        case 1: break;
        case 2: {
            uint16_t *u16 = buf;
            for (int i = 0; i < num; i++) {
                u16[i] = bswap_16(u16[i]);
            }
            break;
        }
        case 4: {
            uint32_t *u32 = buf;
            for (int i = 0; i < num; i++) {
                u32[i] = bswap_32(u32[i]);
            }
            break;
        }
        case 8: {
            uint64_t *u64 = buf;
            for (int i = 0; i < num; i++) {
                u64[i] = bswap_64(u64[i]);
            }
            break;
        }
        default: teal_error("invalid size (%d)", size);
    }
}

void parse_binary(MPI_File file, void *buf, int num, MPI_Datatype type, int swap)
{
    assert((buf || num == 0) && num >= 0);

    if (num == 0) {
        return;
    }

    MPI_Offset offset;
    if (sync.rank == 0) {
        MPI_File_get_position(file, &offset);

        MPI_Status status;
        MPI_File_read_at(file, offset, buf, num, type, &status);

        int count = 0;
        MPI_Get_count(&status, type, &count);
        if (count != num) {
            teal_error("short read (%d / %d)", count, num);
        }

        int size = 0;
        MPI_Type_size(type, &size);

        if (swap) {
            swap_bytes(buf, num, size);
        }

        offset += (MPI_Offset)num * size;
    }

    MPI_Bcast(buf, num, type, 0, sync.comm);

    MPI_Bcast(&offset, 1, MPI_OFFSET, 0, sync.comm);
    MPI_File_seek(file, offset, MPI_SEEK_SET);
}

void parse_binary_split(MPI_File file, void *buf, int num, MPI_Datatype type, int swap)
{
    assert((buf || num == 0) && num >= 0);

    MPI_Offset offset;
    MPI_File_get_position(file, &offset);

    int size = 0;
    MPI_Type_size(type, &size);

    MPI_Offset prefix = num;
    sync_prefix(&prefix, 1, MPI_OFFSET);

    if (num > 0) {
        MPI_Status status;
        MPI_File_read_at(file, offset + (prefix * size), buf, num, type, &status);

        int count = 0;
        MPI_Get_count(&status, type, &count);
        if (count != num) {
            teal_error("short read (%d / %d)", count, num);
        }

        if (swap) {
            swap_bytes(buf, num, size);
        }
    }

    MPI_Offset total = num;
    sync_sum(&total, 1, MPI_OFFSET);

    MPI_File_seek(file, offset + (total * size), MPI_SEEK_SET);
}

void parse(MPI_File file, void *buf, int num, MPI_Datatype type, int mode)
{
    if (mode & BINARY) {
        if (mode & SPLIT) {
            parse_binary_split(file, buf, num, type, mode & SWAP);
            return;
        }
        parse_binary(file, buf, num, type, mode & SWAP);
        return;
    }
    if (mode & SPLIT) {
        parse_ascii_split(file, buf, num, type);
        return;
    }
    parse_ascii(file, buf, num, type);
}
