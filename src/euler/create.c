#include <string.h>

#include "boundary.h"
#include "core/memory.h"
#include "core/utils.h"
#include "euler.h"
#include "flux.h"
#include "timestep.h"
#include "update.h"

static void create_scalar(Equations *eqns, long scalar, const char *name, double value);

static Flux *select_flux(const char *name);

static ApplyBC *select_bc(const char *name);

Equations euler_create(const Mesh *mesh, long space_order)
{
    Equations eqns = {0};

    eqns.n_vars = N_VARS;
    eqns.n_scalars = N_SCALARS;
    strcpy(eqns.name, "Euler");
    eqns.mesh = mesh;

    eqns.vars.name = memory_calloc(N_VARS, sizeof(*eqns.vars.name));
    strcpy(eqns.vars.name[D], "density");
    strcpy(eqns.vars.name[U], "velocity-x");
    strcpy(eqns.vars.name[V], "velocity-y");
    strcpy(eqns.vars.name[W], "velocity-z");
    strcpy(eqns.vars.name[P], "pressure");
    strcpy(eqns.vars.name[DU], "momentum-x");
    strcpy(eqns.vars.name[DV], "momentum-y");
    strcpy(eqns.vars.name[DW], "momentum-z");
    strcpy(eqns.vars.name[DE], "energy");

    eqns.scalar.name = memory_calloc(N_SCALARS, sizeof(*eqns.scalar.name));
    eqns.scalar.value = memory_calloc(N_SCALARS, sizeof(*eqns.scalar.value));
    create_scalar(&eqns, GAMMA, "heat capacity ratio", 1.4);

    eqns.flux.select = select_flux;
    eqns.bc.select = select_bc;

    eqns.timestep = timestep;
    equations_set_flux(&eqns, "hllc");
    eqns.boundary = prim_to_cons;
    eqns.advance = cons_to_prim;

    equations_create(&eqns, space_order);

    return eqns;
}

static void create_scalar(Equations *eqns, long scalar, const char *name, double value)
{
    strcpy(eqns->scalar.name[scalar], name);
    eqns->scalar.value[scalar] = value;
}

static Flux *select_flux(const char *name)
{
    if (!strcmp(name, "godunov"))
        return godunov;
    else if (!strcmp(name, "roe"))
        return roe;
    else if (!strcmp(name, "hll"))
        return hll;
    else if (!strcmp(name, "hllc"))
        return hllc;
    else if (!strcmp(name, "hlle"))
        return hlle;
    else if (!strcmp(name, "lxf"))
        return lxf;
    else
        error("unsupported flux function '%s'", name);
}

static ApplyBC *select_bc(const char *name)
{
    if (!strcmp(name, "symmetry") || !strcmp(name, "wall") || !strcmp(name, "slipwall"))
        return symmetry;
    else if (!strcmp(name, "supersonic inflow"))
        return supersonic_inflow;
    else if (!strcmp(name, "supersonic outflow"))
        return supersonic_outflow;
    else if (!strcmp(name, "subsonic inflow"))
        return subsonic_inflow;
    else if (!strcmp(name, "subsonic outflow"))
        return subsonic_outflow;
    else if (!strcmp(name, "farfield"))
        return farfield;
    else
        error("unsupported boundary condition '%s'", name);
}
