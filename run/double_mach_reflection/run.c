#define _GNU_SOURCE
#include <math.h>

#include "euler.h"

scalar alpha = 30 * M_PI / 180;
Euler upstream = {.density = 1.4, .pressure = 1};
Euler downstream = {.density = 8, .velocity = {.x = 7.1447096, .y = -4.125}, .pressure = 116.5};
scalar shift = 1.0 / 6, Wx = 8.660254, Wy = -5;
Compute dmr;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    vector min_coord = {.x = 0, .y = 0};
    vector max_coord = {.x = 3.25, .y = 1};
    tuple num_cells = {.x = 1040, .y = 320};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, 0);

    vector root = {.x = shift};
    vector normal = {.x = 1};
    mesh_split(mesh, "bottom", root, normal);

    mesh_generate(mesh);
    mesh_summary(mesh);

    Equations *eqns = euler_create(mesh);
    equations_set_boundary_condition(eqns, "left", "supersonic inflow", &downstream, 0);
    equations_set_boundary_condition(eqns, "right", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "bottom-a", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "bottom-b", "slipwall", 0, 0);
    equations_set_boundary_condition(eqns, "top", "custom", 0, dmr);
    equations_set_initial_condition(eqns, "domain", dmr, 0);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 0.2);
    simulation_set_out_time(sim, 0.02);
    simulation_summary(sim);

    simulation_run(sim);

    teal_finalize();
}

void dmr(void *variable_, const scalar *property, vector center, scalar time, const void *ctx_)
{
    Euler *variable = variable_;
    scalar angle = atan2(center.x - shift - (Wx * time), center.y - (Wy * time));
    *variable = (angle >= alpha) ? upstream : downstream;
}
