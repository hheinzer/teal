#include "utils.h"

long x__min_long(long a, long b)
{
    return (a < b ? a : b);
}
double x__min_double(double a, double b)
{
    return (a < b ? a : b);
}

long x__max_long(long a, long b)
{
    return (a > b ? a : b);
}
double x__max_double(double a, double b)
{
    return (a > b ? a : b);
}

long x__sq_long(long a)
{
    return a * a;
}
double x__sq_double(double a)
{
    return a * a;
}
