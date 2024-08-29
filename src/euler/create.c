#include "boundary.h"
#include "euler.h"
#include "flux.h"
#include "timestep.h"
#include "update.h"

Equations *euler_create(const Mesh *mesh)
{
    Equations *eqns = equations_create(mesh, "Euler");

    equations_set_timestep(eqns, euler_timestep);

    equations_set_update(eqns, euler_prim_to_cons, euler_cons_to_prim);

    const String vars[] = {[D] = "density",     [DU] = "momentum-x", [DV] = "momentum-y",
                           [DW] = "momentum-z", [DE] = "energy",     [U] = "velocity-x",
                           [V] = "velocity-y",  [W] = "velocity-z",  [P] = "pressure"};
    equations_create_vars(eqns, vars, N_CONS, N_VARS);

    const String name[] = {[GAMMA] = "heat capacity ratio"};
    equations_create_scalar(eqns, name, N_SCALARS);
    equations_set_scalar(eqns, GAMMA, 1.4);

    equations_set_select_convective_flux(eqns, euler_conv_flux);
    equations_set_convective_flux(eqns, "hllc");

    equations_set_space_order(eqns, 2);
    equations_set_limiter(eqns, "minmod", 0);

    equations_set_select_boundary_condition(eqns, euler_bc);

    return eqns;
}
