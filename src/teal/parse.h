#pragma once

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

typedef struct {
    const char *fname;
    FILE *stream;
} File;

typedef enum { ASCII, BINARY } Mode;

typedef enum { I8, I16, I32, I64, U8, U16, U32, U64, F32, F64, STR } Type;

File parse_open(const char *fname);

/* Return byte offset of file on rank `0`; `-1` on other ranks. */
long parse_get_offset(File file);

/* Seek file to offset on rank `0`; no-op on other ranks. */
void parse_set_offset(File file, long offset);

/* Parse `num` ASCII tokens from file into data. For `type == STR`, `num` is the destination
 * buffer size; result NUL-padded. Data is read on rank `0` and broadcast to other ranks. */
void parse_ascii(Type type, void *data, long num, File file);

/* Parse `num` binary items from file into data; optionally perform a byte swap. Data is read on
 * rank `0` and broadcast to other ranks. */
void parse_binary(Type type, void *data, long num, bool swap, File file);

/* Dispatch `parse_ascii()` or `parse_binary()` based on mode; `swap` is ignored for ASCII. */
void parse(Mode mode, Type type, void *data, long num, bool swap, File file);

/* Parse `num` groups of `len` ASCII tokens with a stride of `stride` tokens into data. Set file
 * position to `beg + len` if `stride > len`, else to end of the full block excluding the last
 * `stride - len` tokens. Return that end offset. Data is read on rank `0` and then sent to the
 * destination rank. */
long parse_split_ascii(Type type, void *data, long num, long len, long stride, File file);

/* Parse `num` groups of `len` binary items with a stride of `stride` bytes into data; optionally
 * perform a byte swap. Set the file position to `beg + len * sizeof(type)` if `stride > len *
 * sizeof(type)`, else to end of the full block excluding the last `stride - len * sizeof(type)`
 * bytes. Return that end offset. Data is read independently by all ranks. */
long parse_split_binary(Type type, void *data, long num, long len, long stride, bool swap,
                        File file);

/* Dispatch `parse_split_ascii()` or `parse_split_binary()` based on mode; `swap` is ignored for
 * ASCII. In ASCII `stride` counts tokens, in BINARY `stride` counts bytes. */
long parse_split(Mode mode, Type type, void *data, long num, long len, long stride, bool swap,
                 File file);

/* Return value at `idx` from data interpreted as type, cast to long (rounded for `F32/F64`). */
long parse_data_to_long(Type type, const void *data, long idx);

void parse_close(File file);
