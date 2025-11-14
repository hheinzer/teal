#pragma once

#include "simulation.h"

typedef struct {
    number num;
    scalar *alpha;
} Lserk;

typedef struct {
    number num;
    scalar *alpha;
} Ssprk;

typedef struct {
    scalar tol_newton;
    scalar tol_krylov;
    number dim_krylov;
} ImplicitEuler;

Advance explicit_euler, lserk, ssprk, implicit_euler;
