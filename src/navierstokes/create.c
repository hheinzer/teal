#include <string.h>

#include "boundary.h"
#include "core/memory.h"
#include "core/utils.h"
#include "euler/boundary.h"
#include "euler/flux.h"
#include "euler/update.h"
#include "flux.h"
#include "navierstokes.h"
#include "timestep.h"

static void create_scalar(Equations *eqns, long scalar, const char *name, double value);

static ConvFlux *select_conv_flux(const char *name);

static ViscFlux *select_visc_flux(const char *name);

static ApplyBC *select_bc(const char *name);

Equations navierstokes_create(const Mesh *mesh, long space_order)
{
    if (space_order < 2) error("unsupported space order '%ld'", space_order);

    Equations eqns = {0};

    eqns.n_cons = N_CONS;
    eqns.n_vars = N_VARS;
    eqns.n_scalars = N_SCALARS;
    strcpy(eqns.name, "Navier-Stokes");
    eqns.mesh = mesh;

    eqns.vars.name = memory_calloc(N_VARS, sizeof(*eqns.vars.name));
    strcpy(eqns.vars.name[D], "density");
    strcpy(eqns.vars.name[DU], "momentum-x");
    strcpy(eqns.vars.name[DV], "momentum-y");
    strcpy(eqns.vars.name[DW], "momentum-z");
    strcpy(eqns.vars.name[DE], "energy");
    strcpy(eqns.vars.name[U], "velocity-x");
    strcpy(eqns.vars.name[V], "velocity-y");
    strcpy(eqns.vars.name[W], "velocity-z");
    strcpy(eqns.vars.name[P], "pressure");

    eqns.scalar.name = memory_calloc(N_SCALARS, sizeof(*eqns.scalar.name));
    eqns.scalar.value = memory_calloc(N_SCALARS, sizeof(*eqns.scalar.value));
    create_scalar(&eqns, GAMMA, "heat capacity ratio", 1.4);
    create_scalar(&eqns, PRANDTL, "Prandtl number", 0.72);
    create_scalar(&eqns, MU, "dynamic viscosity", 0);

    eqns.flux.select_conv = select_conv_flux;
    eqns.flux.select_visc = select_visc_flux;
    eqns.bc.select = select_bc;

    eqns.timestep = timestep_visc;
    equations_set_convective_flux(&eqns, "hllc");
    equations_set_viscous_flux(&eqns, "central");
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

static ConvFlux *select_conv_flux(const char *name)
{
    if (!strcmp(name, "godunov")) return godunov;
    if (!strcmp(name, "roe")) return roe;
    if (!strcmp(name, "hll")) return hll;
    if (!strcmp(name, "hllc")) return hllc;
    if (!strcmp(name, "hlle")) return hlle;
    if (!strcmp(name, "lxf")) return lxf;
    error("unsupported flux function '%s'", name);
}

static ViscFlux *select_visc_flux(const char *name)
{
    if (!strcmp(name, "central")) return central;
    error("unsupported flux function '%s'", name);
}

static ApplyBC *select_bc(const char *name)
{
    if (!strcmp(name, "symmetry") || !strcmp(name, "slipwall")) return symmetry;
    if (!strcmp(name, "adiabatic wall") || !strcmp(name, "wall")) return adiabatic_wall;
    if (!strcmp(name, "supersonic inflow")) return supersonic_inflow;
    if (!strcmp(name, "supersonic outflow")) return supersonic_outflow;
    if (!strcmp(name, "subsonic inflow")) return subsonic_inflow;
    if (!strcmp(name, "subsonic outflow")) return subsonic_outflow;
    if (!strcmp(name, "farfield")) return farfield;
    error("unsupported boundary condition '%s'", name);
}
