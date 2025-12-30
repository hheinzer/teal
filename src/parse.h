#pragma once

#include <mpi.h>

typedef enum {
    ASCII = 1,
    BINARY = 2,
    SPLIT = 4,
    SWAP = 8,
} Mode;

typedef struct Parse {
    MPI_File handle;
    MPI_Offset offset;
} Parse;

// Initialize a parser for the given MPI file.
Parse *parse_init(const char *fname);

// Destroy a parser and close its MPI file.
void parse_deinit(Parse *file);

// Read one string token into a fixed buffer.
int parse_string(Parse *file, char *str, int size);

// Read ASCII tokens into typed elements.
int parse_ascii(Parse *file, void *buf, int num, MPI_Datatype datatype);

// Read ASCII tokens into typed elements with a split count per rank.
int parse_ascii_split(Parse *file, void *buf, int num, MPI_Datatype datatype);

// Read binary data into typed elements.
int parse_binary(Parse *file, void *buf, int num, MPI_Datatype datatype, int mode);

// Read binary data into typed elements with a split count per rank.
int parse_binary_split(Parse *file, void *buf, int num, MPI_Datatype datatype, int mode);

// Dispatch parsing based on the provided mode.
int parse(Parse *file, void *buf, int num, MPI_Datatype datatype, int mode);
