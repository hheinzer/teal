#include "navier_stokes.h"

scalar mach = 0.1, reynolds = 1000;
Compute moving_wall;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    vector min_coord = {.x = 0, .y = 0};
    vector max_coord = {.x = 1, .y = 1};
    tuple num_cells = {.x = 100, .y = 100};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, 0);
    mesh_generate(mesh);
    mesh_summary(mesh);

    NavierStokes state = {.density = 1.4, .pressure = 1};

    Equations *eqns = navier_stokes_create(mesh);
    equations_set_limiter(eqns, "venkatakrishnan", 1);
    equations_set_boundary_condition(eqns, "left", "wall", 0, 0);
    equations_set_boundary_condition(eqns, "right", "wall", 0, 0);
    equations_set_boundary_condition(eqns, "bottom", "wall", 0, 0);
    equations_set_boundary_condition(eqns, "top", "custom", 0, moving_wall);
    equations_set_initial_state(eqns, "domain", &state);
    equations_set_property(eqns, "dynamic viscosity", 1.4 * mach / reynolds);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_iter(sim, 10000);
    simulation_set_out_iter(sim, 100);
    simulation_set_termination(sim, "momentum-x", 1e-5);
    simulation_set_advance(sim, "implicit euler", 100, 0);
    simulation_summary(sim);

    simulation_run(sim);

    teal_finalize();
}

void moving_wall(void *variable_, const scalar *property, vector center, scalar time,
                 const void *ctx_)
{
    NavierStokes *ghost = variable_;
    const NavierStokes *inner = ctx_;
    ghost->density = inner->density;
    ghost->velocity.x = (2 * mach) - inner->velocity.x;
    ghost->velocity.y = -inner->velocity.y;
    ghost->velocity.z = -inner->velocity.z;
    ghost->pressure = inner->pressure;
}
