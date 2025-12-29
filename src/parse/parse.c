#include "parse.h"

#include <assert.h>

#include "teal.h"

int parse(Parse *file, void *buf, int num, MPI_Datatype datatype, int mode)
{
    switch (mode) {
        case ASCII: return parse_ascii(file, buf, num, datatype);
        case ASCII | SPLIT: return parse_ascii_split(file, buf, num, datatype);
        case BINARY: return parse_binary(file, buf, num, datatype);
        case BINARY | SPLIT: return parse_binary_split(file, buf, num, datatype);
        default: teal_error("invalid mode (%d)", mode);
    }
}
