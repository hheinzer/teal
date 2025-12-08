#define _GNU_SOURCE
#include <math.h>

#include "navier_stokes.h"

scalar x = 2, mach = 0.1, reynolds = 1000;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    vector min_coord = {.x = -1, .y = 0};
    vector max_coord = {.x = 5, .y = 5};
    tuple num_cells = {.x = 150, .y = 150};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, 0);

    for (long i = 0; i < mesh->nodes.num; i++) {
        vector *coord = &mesh->nodes.coord[i];
        scalar r = 3;
        coord->y = 5 * expm1(r * coord->y / 5) / expm1(r);
    }

    vector root = {.x = 0, .y = 0};
    vector normal = {.x = 1, .y = 0};
    mesh_split(mesh, "bottom", root, normal);

    mesh_generate(mesh);
    mesh_summary(mesh);

    NavierStokes farfield = {.density = 1.4, .velocity = {.x = mach}, .pressure = 1};

    Equations *eqns = navier_stokes_create(mesh);
    equations_set_limiter(eqns, "venkatakrishnan", 1);
    equations_set_boundary_condition(eqns, "left", "farfield", &farfield, 0);
    equations_set_boundary_condition(eqns, "right", "farfield", &farfield, 0);
    equations_set_boundary_condition(eqns, "bottom-a", "slipwall", 0, 0);
    equations_set_boundary_condition(eqns, "bottom-b", "adiabatic wall", 0, 0);
    equations_set_boundary_condition(eqns, "top", "farfield", &farfield, 0);
    equations_set_initial_state(eqns, "domain", &farfield);
    equations_set_property(eqns, "dynamic viscosity", 1.4 * mach * x / reynolds);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_iter(sim, 1000);
    simulation_set_out_iter(sim, 50);
    simulation_set_termination(sim, "momentum-x", 1e-5);
    simulation_set_advance(sim, "implicit euler", 100, 0);
    simulation_summary(sim);

    simulation_run(sim);

    teal_finalize();
}
