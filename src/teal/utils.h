#pragma once

#define alias(a, b) __typeof__(*b) *a = b

#define min(a, b) _Generic(a, long: x__min_long, double: x__min_double)(a, b)
long x__min_long(long a, long b);
double x__min_double(double a, double b);

#define max(a, b) _Generic(a, long: x__max_long, double: x__max_double)(a, b)
long x__max_long(long a, long b);
double x__max_double(double a, double b);

#define sq(a) _Generic(a, long: x__sq_long, double: x__sq_double)(a)
long x__sq_long(long a);
double x__sq_double(double a);
