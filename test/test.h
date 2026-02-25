#pragma once

#define TEST(expr)                                                                  \
    if (!(expr)) {                                                                  \
        teal2_error("`%s` failed (%s:%d %s)", #expr, __FILE__, __LINE__, __func__); \
    }
