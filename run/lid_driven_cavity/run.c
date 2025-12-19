#include <math.h>

#include "navier_stokes.h"
#include "teal/utils.h"

scalar mach = 0.1, reynolds = 100;
Compute moving_wall;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    vector min_coord = {.x = 0, .y = 0};
    vector max_coord = {.x = 1, .y = 1};
    tuple num_cells = {.x = 100, .y = 100};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, 0);

    for (long i = 0; i < mesh->nodes.num; i++) {
        vector *coord = &mesh->nodes.coord[i];
        scalar beta = 2;
        coord->x = (1 + (tanh(beta * (2 * coord->x - 1)) / tanh(beta))) / 2;
        coord->y = (1 + (tanh(beta * (2 * coord->y - 1)) / tanh(beta))) / 2;
    }

    mesh_generate(mesh);
    mesh_summary(mesh);

    NavierStokes state = {.density = 1, .pressure = 1 / (1.4 * sq(mach))};

    Equations *eqns = navier_stokes_create(mesh);
    equations_set_limiter(eqns, "venkatakrishnan", 1);
    equations_set_boundary_condition(eqns, "left", "wall", 0, 0);
    equations_set_boundary_condition(eqns, "right", "wall", 0, 0);
    equations_set_boundary_condition(eqns, "bottom", "wall", 0, 0);
    equations_set_boundary_condition(eqns, "top", "custom", 0, moving_wall);
    equations_set_initial_state(eqns, "domain", &state);
    equations_set_property(eqns, "dynamic viscosity", 1 / reynolds);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_iter(sim, 10000);
    simulation_set_out_iter(sim, 1000);
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
    ghost->velocity.x = 2 - inner->velocity.x;
    ghost->velocity.y = -inner->velocity.y;
    ghost->velocity.z = -inner->velocity.z;
    ghost->pressure = inner->pressure;
}
