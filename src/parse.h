#pragma once

#include <mpi.h>

enum ParseModes {
    BINARY = 1 << 0,
    SPLIT = 1 << 1,
    SWAP = 1 << 2,
};

// Open a file in read-only mode.
MPI_File parse_open(const char *fname);

// Close a file handle.
void parse_close(MPI_File file);

// Read a string token into `str` (zero-padded to `size`).
int parse_string(MPI_File file, char *str, int size);

// Parse `num` ASCII tokens into `buf` on rank 0 and broadcast.
void parse_ascii(MPI_File file, void *buf, int num, MPI_Datatype type);

// Parse `num` ASCII tokens per rank in rank order.
void parse_ascii_split(MPI_File file, void *buf, int num, MPI_Datatype type);

// Read `num` binary items on rank 0 (optional byte swap) and broadcast.
void parse_binary(MPI_File file, void *buf, int num, MPI_Datatype type, int swap);

// Read `num` binary items per rank (optional byte swap).
void parse_binary_split(MPI_File file, void *buf, int num, MPI_Datatype type, int swap);

// Dispatch to ASCII/BINARY and SPLIT variants based on `mode`.
void parse(MPI_File file, void *buf, int num, MPI_Datatype type, int mode);
