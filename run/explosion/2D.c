#include <math.h>

#include "euler.h"

Euler inner = {.density = 1, .pressure = 1};
Euler outer = {.density = 0.125, .pressure = 0.1};
Compute initial;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    vector min_coord = {.x = -1, .y = -1};
    vector max_coord = {.x = 1, .y = 1};
    tuple num_cells = {.x = 100, .y = 100};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, 0);
    mesh_generate(mesh);
    mesh_summary(mesh);

    Equations *eqns = euler_create(mesh);
    equations_set_boundary_condition(eqns, "left", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "right", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "bottom", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "top", "supersonic outflow", 0, 0);
    equations_set_initial_condition(eqns, "domain", initial, 0);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 0.25);
    simulation_summary(sim);

    simulation_run(sim);

    teal_finalize();
}

void initial(void *variable_, const scalar *property, vector center, scalar time, const void *ctx_)
{
    Euler *variable = variable_;
    scalar radius = hypot(center.x, center.y);
    *variable = (radius <= 0.4) ? inner : outer;
}
