#include <assert.h>

#include "euler.h"

Equations *euler_create(const Mesh *mesh)
{
    assert(mesh);

    Equations *eqns = equations_create(mesh, "euler");

    long size[5] = {1, 3, 1, 3, 1};
    const char *variable[5] = {"density", "momentum", "energy", "velocity", "pressure"};
    equations_create_variables(eqns, size, variable, euler_conserved, euler_primitive, 3, 5);

    const char *property[1] = {"heat capacity ratio"};
    const scalar gamma = 1.4;  // air
    equations_create_properties(eqns, property, &gamma, 1);

    equations_set_space_order(eqns, 2);
    equations_set_timestep(eqns, euler_timestep);
    equations_set_boundary_select(eqns, euler_boundary);
    equations_set_convective_select(eqns, euler_convective);
    equations_set_limiter(eqns, "minmod", 0);
    equations_set_convective_flux(eqns, "hllc");

    return eqns;
}
