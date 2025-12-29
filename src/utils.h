#pragma once

#define trap(rnk)                                 \
    do {                                          \
        if (sync.rank == (rnk)) __builtin_trap(); \
    } while (0)

#define breakpoint(rnk)                                       \
    do                                                        \
        if (sync.rank == (rnk)) __asm__ __volatile__("int3"); \
    }                                                         \
    while (0)

#define sq(val) ((val) * (val))
#define cb(val) ((val) * (val) * (val))

#define min(lhs, rhs) (((lhs) < (rhs)) ? (lhs) : (rhs))
#define max(lhs, rhs) (((lhs) > (rhs)) ? (lhs) : (rhs))

// Compare floating point values with tolerances.
int isclose(double lhs, double rhs);
