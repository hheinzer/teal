#include <assert.h>
#include <limits.h>
#include <stddef.h>
#include <string.h>

#include "parse.h"
#include "space.h"
#include "sync.h"
#include "teal.h"
#include "utils.h"

// Return non zero if the character is a quote delimiter.
static int isquote(char chr)
{
    switch (chr) {
        case '"':
        case '\'':
        case '`': return 1;
        default: return 0;
    }
}

// Find the matching quote or return end if not found.
static char *find_quote(char *beg, const char *end)
{
    char quote = *beg;
    beg += 1;
    while (beg < end && *beg != quote) {
        beg += 1;
    }
    return beg;
}

// Parse one token from a buffer and advance the file offset.
static int read_token(Parse *file, char *str, int size, int count)
{
    // locate the token bounds in the chunk
    const char *end = str + count;
    char *beg = skip_space(str, end);
    if (beg >= end) {
        teal_error("invalid string (only spaces)");
    }

    if (isquote(*beg)) {
        // quoted tokens stop at the matching quote
        const char *quote = find_quote(beg, end);
        if (quote >= end) {
            teal_error("unterminated quoted string");
        }
        end = quote;
        beg += 1;
    }
    else {
        // unquoted tokens stop at the next delimiter
        const char *space = find_space(beg, end);
        if (space >= end && count == size) {
            teal_error("missing delimiter (buffer too small)");
        }
        end = space;
    }

    // advance past the token and past a delimiter if present
    file->offset += end - str;
    if (end < str + count) {
        file->offset += 1;
    }

    // pack the token into the destination buffer
    ptrdiff_t diff = end - beg;
    assert(inrange(0, diff, size));
    int len = (int)diff;
    memmove(str, beg, (size_t)len);
    memset(str + len, 0, (size_t)(size - len));
    return len;
}

int parse_string(Parse *file, char *str, int size)
{
    assert(file && (str || size == 0) && size >= 0);
    if (size == 0) {
        return 0;
    }
    int len = 0;
    if (sync.rank == 0) {
        // read a fixed chunk for one token
        MPI_Status status;
        MPI_File_read_at(file->handle, file->offset, str, size, MPI_CHAR, &status);
        int count = 0;
        MPI_Get_count(&status, MPI_CHAR, &count);
        if (count <= 0) {
            memset(str, 0, (size_t)size);
        }
        else {
            len = read_token(file, str, size, count);
        }
    }
    MPI_Bcast(&len, 1, MPI_INT, 0, sync.comm);
    MPI_Bcast(str, size, MPI_CHAR, 0, sync.comm);
    MPI_Bcast(&file->offset, 1, MPI_OFFSET, 0, sync.comm);
    return len;
}
