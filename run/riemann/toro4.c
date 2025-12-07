#include <stdio.h>
#include <stdlib.h>

#include "euler.h"

Euler left = {.density = 5.99924, .velocity = {.x = 19.5975}, .pressure = 460.894};
Euler right = {.density = 5.99242, .velocity = {.x = -6.19633}, .pressure = 46.0950};
scalar x0 = 0.4, max_time = 0.035;
Compute exact;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);
    const char *flux = (argc > 1) ? argv[1] : "hllc";
    long space_order = (argc > 2) ? strtol(argv[2], 0, 10) : 2;

    vector min_coord = {.x = 0};
    vector max_coord = {.x = 1};
    tuple num_cells = {.x = 100};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, 0);
    mesh_generate(mesh);
    mesh_summary(mesh);

    Equations *eqns = euler_create(mesh);
    equations_create_exact_solution(eqns, exact);
    equations_set_convective_flux(eqns, flux);
    equations_set_space_order(eqns, space_order);
    equations_set_boundary_condition(eqns, "left", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "right", "supersonic outflow", 0, 0);
    equations_set_initial_condition(eqns, "domain", exact, 0);
    equations_summary(eqns);

    char prefix[128];
    sprintf(prefix, "%s_%s_%ld", argv[0], flux, space_order);

    Simulation *sim = simulation_create(eqns, prefix);
    simulation_set_max_time(sim, max_time);
    simulation_summary(sim);

    scalar time = simulation_run(sim);
    simulation_error(sim, time, 0);

    teal_finalize();
}

void exact(void *variable_, const scalar *property, vector center, scalar time, const void *ctx_)
{
    Euler *variable = variable_;
    scalar gamma = property[EULER_HEAT_CAPACITY_RATIO];
    scalar location = (center.x - x0) / time;
    *variable = euler_riemann(&left, &right, gamma, location);
}
