#include <ctype.h>

// Skip leading whitespace in a bounded buffer.
static char *skip_space(char *beg, const char *end)
{
    while (beg < end && isspace((unsigned char)*beg)) {
        beg += 1;
    }
    return beg;
}

// Find the next whitespace or return end.
static char *find_space(char *beg, const char *end)
{
    while (beg < end && !isspace((unsigned char)*beg)) {
        beg += 1;
    }
    return beg;
}
