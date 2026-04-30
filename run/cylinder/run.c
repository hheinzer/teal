#include <math.h>

#include "navierstokes.h"
#include "utils.h"

static const double mach = 0.1, reynolds = 250, velocity_r = 0;

static void rotating_wall(void *outer_, const double *property, Vector center, double time,
                          const void *context)
{
    NavierStokesPrimitive *outer = outer_;
    const NavierStokesPrimitive *inner = context;
    double theta = atan2(inner->velocity.y, inner->velocity.x) - (M_PI / 2);
    outer->density = inner->density;
    outer->velocity.x = (2 * velocity_r * cos(theta)) - inner->velocity.x;
    outer->velocity.y = (2 * velocity_r * sin(theta)) - inner->velocity.y;
    outer->velocity.z = -inner->velocity.z;
    outer->pressure = inner->pressure;
}

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);

    Mesh *mesh = mesh_read("run/cylinder/mesh.msh");
    mesh_generate(mesh);
    mesh_summary(mesh);

    NavierStokesPrimitive farfield = {
        .density = 1,
        .velocity = {.x = 1},
        .pressure = 1 / (1.4 * sq(mach)),
    };

    Equations *eqns = navierstokes_create(mesh);
    equations_set_limiter(eqns, "venkatakrishnan", 1);
    equations_set_boundary_condition(eqns, "wall", "custom", rotating_wall, 0);
    equations_set_boundary_condition(eqns, "farfield", "farfield", 0, &farfield);
    equations_set_initial_condition(eqns, "domain", 0, &farfield);
    equations_set_property(eqns, "dynamic viscosity", 1 / reynolds);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 100);
    simulation_set_out_time(sim, 1);
    simulation_summary(sim);

    simulation_run(sim);

    simulation_destroy(sim);
    equations_destroy(eqns);
    mesh_destroy(mesh);

    teal_deinit();
}
