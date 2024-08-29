#include "euler.h"

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    Mesh *mesh = mesh_read("run/ffs/quads.geo");
    mesh_generate(&mesh);
    mesh_print(mesh);

    const double state[N_VARS] = {[D] = 1.4, [U] = 3, [P] = 1};
    Equations *eqns = euler_create(mesh);
    equations_set_boundary_condition(eqns, "inlet", "supersonic inflow", state, 0);
    equations_set_boundary_condition(eqns, "outlet", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "symmetry", "symmetry", 0, 0);
    equations_set_boundary_condition(eqns, "wall", "slipwall", 0, 0);
    equations_set_initial_state(eqns, state);
    equations_print(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 2);
    simulation_print(sim);
    simulation_run(sim);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);

    teal_finalize();
}
