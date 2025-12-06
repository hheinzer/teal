#include <string.h>

#include "euler.h"
#include "navier_stokes.h"

static void adiabatic_wall(void *ghost_, const void *inner_, const void *reference_,
                           const scalar *property, const Basis *basis)
{
    (void)reference_;
    (void)property;
    (void)basis;
    NavierStokes *ghost = ghost_;
    const NavierStokes *inner = inner_;

    ghost->density = inner->density;
    ghost->velocity.x = -inner->velocity.x;
    ghost->velocity.y = -inner->velocity.y;
    ghost->velocity.z = -inner->velocity.z;
    ghost->pressure = inner->pressure;
}

Boundary *navier_stokes_boundary(const char *name)
{
    if (!strcmp(name, "adiabatic wall") || !strcmp(name, "wall")) {
        return adiabatic_wall;
    }
    return euler_boundary(name);
}
