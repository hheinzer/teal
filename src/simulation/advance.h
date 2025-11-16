#pragma once

#include "simulation.h"

typedef struct {
    number time_order;
    number num_stages;
} Lserk;

Advance euler, midpoint, heun, ralston, ssprk2, ssprk3, rk3, rk4, lserk;
