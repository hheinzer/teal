#include "euler.h"

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);

    Mesh *mesh = mesh_read("run/scramjet_inlet/mesh.msh");
    mesh_generate(mesh);
    mesh_summary(mesh);

    EulerPrimitive inlet = {.density = 1.4, .velocity = {.x = 3}, .pressure = 1};

    Equations *eqns = euler_create(mesh);
    equations_set_boundary_condition(eqns, "inlet", "supersonic inflow", 0, &inlet);
    equations_set_boundary_condition(eqns, "outlet", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "wall", "slipwall", 0, 0);
    equations_set_initial_condition(eqns, "domain", 0, &inlet);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 5);
    simulation_set_out_time(sim, 0.5);
    simulation_summary(sim);

    simulation_run(sim);

    simulation_destroy(sim);
    equations_destroy(eqns);
    mesh_destroy(mesh);

    teal_deinit();
}
