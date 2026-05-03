#include "navierstokes.h"

static const double mach = 0.3, reynolds = 1000;

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);

    Mesh *mesh = mesh_read("run/sphere/mesh.msh");
    mesh_generate(mesh);
    mesh_summary(mesh);

    NavierStokesPrimitive farfield = {
        .density = 1, .velocity = {.x = 1}, .pressure = 1 / (1.4 * sq(mach))};

    Equations *eqns = navierstokes_create(mesh);
    equations_set_limiter(eqns, "venkatakrishnan", 1);
    equations_set_boundary_condition(eqns, "wall", "wall", 0, 0);
    equations_set_boundary_condition(eqns, "farfield", "farfield", 0, &farfield);
    equations_set_initial_condition(eqns, "domain", 0, &farfield);
    equations_set_property(eqns, "dynamic viscosity", 1 / reynolds);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 100);
    simulation_set_out_time(sim, 0.1);
    simulation_summary(sim);

    simulation_run(sim);

    simulation_destroy(sim);
    equations_destroy(eqns);
    mesh_destroy(mesh);

    teal_deinit();
}
