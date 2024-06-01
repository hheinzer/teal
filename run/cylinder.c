#include "teal.h"

int main(int argc, char **argv) {
    teal_init(argc, argv);

    const double state[] = {1.4, 0.1, 0, 1};
    Mesh mesh = mesh_read("geo/cylinder.geo");
    mesh_set_boundary_condition(&mesh, "wall", "wall", 0, 0);
    mesh_set_boundary_condition(&mesh, "farfield", "characteristic", state, 0);
    mesh_generate(&mesh);
    mesh_print(&mesh);

    Equations eqns = euler_create(&mesh, 0);
    equations_set_initial_state(&eqns, 4, (long[]){D, U, V, P}, state);
    euler_compute_conserved(&eqns);
    euler_print(&eqns);

    Simulation sim = simulation_create(argv[0], &eqns);
    sim.abort_residual = 1e-4;
    sim.output_iter = 1000;
    simulation_print(&sim);
    simulation_run(&sim);
    airfoil_polar(&sim, "wall", P, state, 2);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);
    teal_finalize();
}
