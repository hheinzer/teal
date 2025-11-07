#include "assert.h"

#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include "sync.h"

#define UNW_LOCAL_ONLY
#include <libunwind.h>

static bool append(char *buf, int size, int *pos, const char *fmt, ...)
{
    if (*pos >= size) {
        buf[size - 1] = 0;
        return false;
    }

    va_list args;
    va_start(args, fmt);
    int rem = size - *pos;
    int num = vsnprintf(&buf[*pos], rem, fmt, args);
    va_end(args);

    if (num < 0) {
        return false;
    }

    if (num >= rem) {
        *pos = size - 1;
        buf[*pos] = 0;
        return false;
    }

    *pos += num;
    return true;
}

void x__assert_fail(const char *file, int line, const char *func, const char *expr)
{
    char buf[4096];
    int pos = 0;
    if (!append(buf, sizeof(buf), &pos, "[%d] %s:%d: %s: Assertion `%s` failed.\n", sync.rank, file,
                line, func, expr)) {
        goto out;
    }

    unw_context_t ctx;
    unw_cursor_t cur;
    if (!unw_getcontext(&ctx) && !unw_init_local(&cur, &ctx)) {
        int frame = 0;
        while (unw_step(&cur) > 0) {
            char name[128];
            unw_word_t off = 0;
            if (unw_get_proc_name(&cur, name, sizeof(name), &off)) {
                strcpy(name, "???");
                off = 0;
            }
            unw_word_t iptr = 0;
            unw_get_reg(&cur, UNW_REG_IP, &iptr);
            if (!append(buf, sizeof(buf), &pos, "\t %2d. %-30s (+0x%lx) [0x%lx]\n", frame++, name,
                        off, iptr)) {
                break;
            }
        }
    }

out:
    if (pos > 0) {
        fwrite(buf, sizeof(*buf), pos, stderr);
        fflush(stderr);
    }
    sync_abort();
}
