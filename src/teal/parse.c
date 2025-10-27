#include "parse.h"

#include <byteswap.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "arena.h"
#include "assert.h"
#include "sync.h"

enum { MAX_TOKEN = 4096 };

static const MPI_Datatype datatype_of[] = {
    [I8] = MPI_INT8_T,  [I16] = MPI_INT16_T,  [I32] = MPI_INT32_T,  [I64] = MPI_INT64_T,
    [U8] = MPI_UINT8_T, [U16] = MPI_UINT16_T, [U32] = MPI_UINT32_T, [U64] = MPI_UINT64_T,
    [F32] = MPI_FLOAT,  [F64] = MPI_DOUBLE,   [STR] = MPI_CHAR,
};

static const long size_of[] = {
    [I8] = sizeof(int8_t),    [I16] = sizeof(int16_t),  [I32] = sizeof(int32_t),
    [I64] = sizeof(int64_t),  [U8] = sizeof(uint8_t),   [U16] = sizeof(uint16_t),
    [U32] = sizeof(uint32_t), [U64] = sizeof(uint64_t), [F32] = sizeof(float),
    [F64] = sizeof(double),
};

File parse_open(const char *fname)
{
    assert(fname);
    File file = {0};
    if (sync.rank == 0) {
        file.stream = fopen(fname, "rb");
        assert(file.stream);
    }
    MPI_File_open(sync.comm, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &file.handle);
    return file;
}

long parse_get_offset(File file)
{
    long offset = -1;
    if (sync.rank == 0) {
        assert(file.stream);
        offset = ftell(file.stream);
        assert(offset != -1);
    }
    return offset;
}

void parse_set_offset(File file, long offset)
{
    if (sync.rank == 0) {
        assert(file.stream && offset >= 0);
        int ret = fseek(file.stream, offset, SEEK_SET);
        assert(ret == 0);
    }
}

/* Read next ASCII token from stream; unquote string it the `type == STR` and quotes are used. */
static long next(Type type, char *token, long size, FILE *stream)
{
    int chr;
    do {  // skip leading white space
        chr = fgetc(stream);
        assert(chr != EOF);
    } while (isspace((unsigned char)chr));

    long len = 0;
    if (type == STR && (chr == '\"' || chr == '\'')) {
        int quote = chr;
        while (true) {
            chr = fgetc(stream);
            assert(chr != EOF);
            if (chr == quote) {
                break;
            }
            assert(len + 1 < size);
            token[len++] = chr;
        }
        token[len] = 0;
        return len;
    }

    do {
        assert(len + 1 < size);
        token[len++] = chr;
        chr = fgetc(stream);
    } while (chr != EOF && !isspace((unsigned char)chr));

    token[len] = 0;
    return len;
}

/* Convert token string into to the data type at index `idx`. */
static void token_to_data(Type type, void *data, long idx, const char *token)
{
    errno = 0;
    char *end = 0;
    switch (type) {
        case I8: ((int8_t *)data)[idx] = strtol(token, &end, 10); break;
        case I16: ((int16_t *)data)[idx] = strtol(token, &end, 10); break;
        case I32: ((int32_t *)data)[idx] = strtol(token, &end, 10); break;
        case I64: ((int64_t *)data)[idx] = strtoll(token, &end, 10); break;
        case U8: ((uint8_t *)data)[idx] = strtoul(token, &end, 10); break;
        case U16: ((uint16_t *)data)[idx] = strtoul(token, &end, 10); break;
        case U32: ((uint32_t *)data)[idx] = strtoul(token, &end, 10); break;
        case U64: ((uint64_t *)data)[idx] = strtoull(token, &end, 10); break;
        case F32: ((float *)data)[idx] = strtof(token, &end); break;
        case F64: ((double *)data)[idx] = strtod(token, &end); break;
        default: assert(false);
    }
    assert(errno == 0 && end && end != token && *end == 0);
}

void parse_ascii(Type type, void *data, long num, File file)
{
    assert(data && num >= 0);
    if (sync.rank == 0) {
        assert(file.stream);
        char token[MAX_TOKEN];
        if (type == STR) {
            long len = next(type, token, sizeof(token), file.stream);
            assert(len < num);
            char *str = data;
            memcpy(str, token, len);
            memset(str + len, 0, num - len);
        }
        else {
            for (long i = 0; i < num; i++) {
                next(type, token, sizeof(token), file.stream);
                token_to_data(type, data, i, token);
            }
        }
    }
    MPI_Bcast(data, num, datatype_of[type], 0, sync.comm);
}

/* Perform byte swap for data based on the item size. */
static void swap_data(void *data, long num, long size)
{
    switch (size) {
        case 1: break;
        case 2: {
            uint16_t *u16 = data;
            for (long i = 0; i < num; i++) {
                u16[i] = bswap_16(u16[i]);
            }
            break;
        }
        case 4: {
            uint32_t *u32 = data;
            for (long i = 0; i < num; i++) {
                u32[i] = bswap_32(u32[i]);
            }
            break;
        }
        case 8: {
            uint64_t *u64 = data;
            for (long i = 0; i < num; i++) {
                u64[i] = bswap_64(u64[i]);
            }
            break;
        }
        default: assert(false);
    }
}

void parse_binary(Type type, void *data, long num, bool swap, File file)
{
    assert(type != STR && data && num >= 0);
    if (sync.rank == 0) {
        assert(file.stream);
        long read = fread(data, size_of[type], num, file.stream);
        assert(read == num);
        if (swap) {
            swap_data(data, num, size_of[type]);
        }
    }
    MPI_Bcast(data, num, datatype_of[type], 0, sync.comm);
}

void parse(Mode mode, Type type, void *data, long num, bool swap, File file)
{
    switch (mode) {
        case ASCII: parse_ascii(type, data, num, file); break;
        case BINARY: parse_binary(type, data, num, swap, file); break;
        default: assert(false);
    }
}

long parse_split_ascii(Type type, void *data, long num, long len, long stride, File file)
{
    assert(type != STR && data && num >= 0 && len >= 0 && stride >= len);

    long tot = sync_lsum(num);

    int tag_num = sync_tag();
    int tag_data = sync_tag();
    if (sync.rank != 0) {
        MPI_Send(&num, 1, MPI_LONG, 0, tag_num, sync.comm);
        MPI_Recv(data, num * len, datatype_of[type], 0, tag_data, sync.comm, MPI_STATUS_IGNORE);
    }

    long end;
    if (sync.rank == 0) {
        long gap = stride - len;
        long pos = -1;

        char token[MAX_TOKEN];
        long cnt = 0;
        for (long rank = 0; rank < sync.size; rank++) {
            Arena save = arena_save();

            long num_rank = num;
            void *data_rank = data;
            if (rank != 0) {
                MPI_Recv(&num_rank, 1, MPI_LONG, rank, tag_num, sync.comm, MPI_STATUS_IGNORE);
                data_rank = arena_malloc(num_rank * len, size_of[type]);
            }

            for (long i = 0; i < num_rank; i++) {
                for (long j = 0; j < len; j++) {
                    next(type, token, sizeof(token), file.stream);
                    token_to_data(type, data_rank, (i * len) + j, token);
                }
                if (pos < 0) {
                    pos = ftell(file.stream);
                }
                if (++cnt < tot) {  // only skip internal gaps
                    for (long j = 0; j < gap; j++) {
                        next(type, token, sizeof(token), file.stream);
                    }
                }
            }

            if (rank != 0) {
                MPI_Send(data_rank, num_rank * len, datatype_of[type], rank, tag_data, sync.comm);
            }

            arena_load(save);
        }

        end = parse_get_offset(file);

        if (tot > 0 && gap > 0) {
            assert(pos >= 0);
            parse_set_offset(file, pos);
        }
    }

    MPI_Bcast(&end, 1, MPI_LONG, 0, sync.comm);
    return end;
}

long parse_split_binary(Type type, void *data, long num, long len, long stride, bool swap,
                        File file)
{
    assert(type != STR && data && num >= 0 && len >= 0 && stride >= len * size_of[type]);

    long beg = parse_get_offset(file);
    MPI_Bcast(&beg, 1, MPI_LONG, 0, sync.comm);

    long gap = stride - (len * size_of[type]);
    MPI_Offset disp = beg + (sync_lexsum(num) * stride);
    if (gap == 0) {
        MPI_File_read_at(file.handle, disp, data, num * len, datatype_of[type], MPI_STATUS_IGNORE);
    }
    else {
        Arena save = arena_save();

        char (*src)[stride] = arena_malloc(num, sizeof(*src));
        MPI_File_read_at(file.handle, disp, src, num * stride, MPI_BYTE, MPI_STATUS_IGNORE);

        char (*dst)[len * size_of[type]] = data;
        for (long i = 0; i < num; i++) {
            memcpy(&dst[i], &src[i], sizeof(*dst));
        }

        arena_load(save);
    }

    if (swap) {
        swap_data(data, num * len, size_of[type]);
    }

    long tot = sync_lsum(num);
    long end = (tot == 0) ? beg : beg + (tot * stride) - gap;

    parse_set_offset(file, (tot > 0 && gap > 0) ? beg + (len * size_of[type]) : end);
    return end;
}

long parse_split(Mode mode, Type type, void *data, long num, long len, long stride, bool swap,
                 File file)
{
    switch (mode) {
        case ASCII: return parse_split_ascii(type, data, num, len, stride, file);
        case BINARY: return parse_split_binary(type, data, num, len, stride, swap, file);
        default: assert(false);
    }
}

long parse_data_to_long(Type type, const void *data, long idx)
{
    assert(type != STR && data && idx >= 0);
    switch (type) {
        case I8: return ((int8_t *)data)[idx];
        case I16: return ((int16_t *)data)[idx];
        case I32: return ((int32_t *)data)[idx];
        case I64: return ((int64_t *)data)[idx];
        case U8: return ((uint8_t *)data)[idx];
        case U16: return ((uint16_t *)data)[idx];
        case U32: return ((uint32_t *)data)[idx];
        case U64: return ((uint64_t *)data)[idx];
        case F32: return lrintf(((float *)data)[idx]);
        case F64: return lrint(((double *)data)[idx]);
        default: assert(false);
    }
}

void parse_close(File file)
{
    if (sync.rank == 0) {
        assert(file.stream);
        fclose(file.stream);
    }
    MPI_File_close(&file.handle);
}
