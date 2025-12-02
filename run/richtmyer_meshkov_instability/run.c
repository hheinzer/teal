#define _GNU_SOURCE
#include <math.h>

#include "euler.h"

Compute rmi;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    vector min_coord = {0, 0, 0};
    vector max_coord = {60, 15, 15};
    tuple num_cells = {400, 100, 1};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, 0);
    mesh_generate(mesh);
    mesh_summary(mesh);

    Equations *eqns = euler_create(mesh);
    equations_set_limiter(eqns, "venkatakrishnan", 1);
    equations_set_boundary_condition(eqns, "left", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "right", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "bottom", "symmetry", 0, 0);
    equations_set_boundary_condition(eqns, "top", "symmetry", 0, 0);
    equations_set_initial_condition(eqns, "domain", rmi, 0);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 200);
    simulation_set_out_time(sim, 20);
    simulation_summary(sim);

    simulation_run(sim);

    teal_finalize();
}

void rmi(void *variable_, const scalar *property, vector center, scalar time, const void *ctx_)
{
    Euler *variable = variable_;
    variable->pressure = 1;
    if (center.x >= 18 + (2 * cos(2 * M_PI * center.y / 15))) {
        variable->density = 0.25;
    }
    else {
        variable->density = 1;
    }
    if (2 <= center.x && center.x <= 6) {
        variable->density = 4.22;
        variable->pressure = 4.9;
    }
}
