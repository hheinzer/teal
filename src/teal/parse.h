#pragma once

#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>

typedef struct {
    FILE *stream;
    MPI_File handle;
} ParseFile;

typedef enum { ASCII, BINARY } ParseMode;

typedef enum { I8, I16, I32, I64, U8, U16, U32, U64, F32, F64, STR } ParseType;

ParseFile parse_open(const char *fname);
void parse_close(ParseFile file);

long parse_get_offset(ParseFile file);
void parse_set_offset(ParseFile file, long offset);

// Parse `num` ASCII tokens from file into data. For `type == STR`, `num` is the destination
// buffer size; result NUL-padded. Data is read on rank `0` and broadcast to other ranks.
void parse_ascii(ParseType type, void *data, long num, ParseFile file);

// Parse `num` binary items from file into data; optionally perform a byte swap. Data is read on
// rank `0` and broadcast to other ranks.
void parse_binary(ParseType type, void *data, long num, bool swap, ParseFile file);

void parse(ParseMode mode, ParseType type, void *data, long num, bool swap, ParseFile file);

// Parse `num` groups of `len` ASCII tokens with a stride of `stride` tokens into data. Set file
// position to `beg + len` if `stride > len`, else to end of the full block excluding the last
// `stride - len` tokens. Return that end offset. Data is read on rank `0`, split, and each rank
// receives its own block.
long parse_split_ascii(ParseType type, void *data, long num, long len, long stride, ParseFile file);

// Parse `num` groups of `len` binary items with a stride of `stride` bytes into data; optionally
// perform a byte swap. Set the file position to `beg + len * sizeof(type)` if `stride > len *
// sizeof(type)`, else to end of the full block excluding the last `stride - len * sizeof(type)`
// bytes. Return that end offset. Data is read independently by all ranks.
long parse_split_binary(ParseType type, void *data, long num, long len, long stride, bool swap,
                        ParseFile file);

long parse_split(ParseMode mode, ParseType type, void *data, long num, long len, long stride,
                 bool swap, ParseFile file);

long parse_data_to_long(ParseType type, const void *data, long idx);
