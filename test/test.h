#pragma once

#define test(expr) \
    ((expr) ? (void)0 : teal_error("`%s` failed (%s:%d: %s)", #expr, __FILE__, __LINE__, __func__))
