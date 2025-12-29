#include "parse.h"

#include <assert.h>

#include "teal.h"

int parse(Parse *file, void *buf, int num, MPI_Datatype datatype, int mode)
{
    int split = mode & SPLIT;
    int swap = mode & SWAP;
    mode &= ~(SPLIT | SWAP);
    switch (mode) {
        case ASCII:
            return split ? parse_ascii_split(file, buf, num, datatype)
                         : parse_ascii(file, buf, num, datatype);
        case BINARY:
            return split ? parse_binary_split(file, buf, num, datatype, swap)
                         : parse_binary(file, buf, num, datatype, swap);
        default: teal_error("invalid mode (%d)", mode);
    }
}
