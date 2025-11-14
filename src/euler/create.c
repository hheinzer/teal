#include "euler.h"
#include "teal/assert.h"

Equations *euler_create(const Mesh *mesh)
{
    assert(mesh);

    Equations *eqns = equations_create(mesh, "euler");

    const char *variable[5] = {"density", "momentum", "energy", "velocity", "pressure"};
    Type type[5] = {SCALAR, VECTOR, SCALAR, VECTOR, SCALAR};
    equations_create_variables(eqns, variable, type, euler_conserved, euler_primitive, 3, 5);

    const char *property = "heat capacity ratio";
    const scalar gamma = 1.4;  // air
    equations_create_properties(eqns, &property, &gamma, 1);

    equations_set_space_order(eqns, 2);
    equations_set_timestep(eqns, euler_timestep);
    equations_set_boundary_select(eqns, euler_boundary);
    equations_set_convective_select(eqns, euler_convective);
    equations_set_limiter(eqns, "minmod", 0);
    equations_set_convective_flux(eqns, "hllc");

    return eqns;
}
