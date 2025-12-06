#include <assert.h>

#include "euler.h"
#include "navier_stokes.h"

Equations *navier_stokes_create(const Mesh *mesh)
{
    assert(mesh);

    Equations *eqns = equations_create(mesh, "navier-stokes");

    long dim[5] = {1, 3, 1, 3, 1};
    const char *variable[5] = {"density", "momentum", "energy", "velocity", "pressure"};
    equations_create_variables(eqns, dim, variable, euler_conserved, euler_primitive, 3, 5);

    const char *property[3] = {"heat capacity ratio", "dynamic viscosity", "prandtl number"};
    const scalar value[3] = {1.4, 0, 0.71};  // air (except viscosity)
    equations_create_properties(eqns, property, value, 3);

    equations_set_space_order(eqns, 2);
    equations_set_time_step(eqns, navier_stokes_time_step);
    equations_set_boundary_select(eqns, navier_stokes_boundary);
    equations_set_convective_select(eqns, euler_convective);
    equations_set_viscous_select(eqns, navier_stokes_viscous);
    equations_set_limiter(eqns, "minmod", 0);
    equations_set_convective_flux(eqns, "hllc");
    equations_set_viscous_flux(eqns, "central");

    return eqns;
}
