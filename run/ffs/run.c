#include "euler.h"

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    Mesh *mesh = mesh_read("run/ffs/quad.msh");
    mesh_generate(mesh);
    mesh_summary(mesh);

    Euler state = {.density = 1.4, .velocity = {.x = 3}, .pressure = 1};

    Equations *eqns = euler_create(mesh);
    equations_set_boundary_condition(eqns, "inlet", "supersonic inflow", &state, 0);
    equations_set_boundary_condition(eqns, "outlet", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "symmetry", "symmetry", 0, 0);
    equations_set_boundary_condition(eqns, "wall", "slipwall", 0, 0);
    equations_set_initial_state(eqns, "domain", &state);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 2);
    simulation_summary(sim);

    simulation_run(sim);

    teal_finalize();
}
