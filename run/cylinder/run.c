#define _GNU_SOURCE
#include <math.h>

#include "navier_stokes.h"
#include "teal/utils.h"

scalar mach = 0.1, reynolds = 250, velocity_r = 0;
Compute rotating_wall;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    Mesh *mesh = mesh_read("run/cylinder/mesh.msh");
    mesh_generate(mesh);
    mesh_summary(mesh);

    NavierStokes farfield = {.density = 1, .velocity = {.x = 1}, .pressure = 1 / (1.4 * sq(mach))};

    Equations *eqns = navier_stokes_create(mesh);
    equations_set_limiter(eqns, "venkatakrishnan", 1);
    equations_set_boundary_condition(eqns, "wall", "custom", 0, rotating_wall);
    equations_set_boundary_condition(eqns, "farfield", "farfield", &farfield, 0);
    equations_set_initial_state(eqns, "domain", &farfield);
    equations_set_property(eqns, "dynamic viscosity", 1 / reynolds);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 100);
    simulation_set_out_time(sim, 1);
    simulation_summary(sim);

    simulation_run(sim);

    teal_finalize();
}

void rotating_wall(void *variable_, const scalar *property, vector center, scalar time,
                   const void *ctx_)
{
    NavierStokes *ghost = variable_;
    const NavierStokes *inner = ctx_;
    scalar theta = atan2(inner->velocity.y, inner->velocity.x) - M_PI_2;
    ghost->density = inner->density;
    ghost->velocity.x = (2 * velocity_r * cos(theta)) - inner->velocity.x;
    ghost->velocity.y = (2 * velocity_r * sin(theta)) - inner->velocity.y;
    ghost->velocity.z = -inner->velocity.z;
    ghost->pressure = inner->pressure;
}
