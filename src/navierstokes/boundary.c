#include <string.h>

#include "euler.h"
#include "navierstokes.h"

static void adiabatic_wall(void *outer_, const void *inner_, const double *property, Matrix basis,
                           const void *context)
{
    (void)property;
    (void)context;
    (void)basis;
    NavierStokesPrimitive *outer = outer_;
    const NavierStokesPrimitive *inner = inner_;

    outer->density = inner->density;
    outer->velocity = vector_mul(-1, inner->velocity);
    outer->pressure = inner->pressure;
}

Boundary *navierstokes_boundary(const char *name)
{
    if (!strcmp(name, "adiabatic wall") || !strcmp(name, "wall")) {
        return adiabatic_wall;
    }
    return euler_boundary(name);
}
