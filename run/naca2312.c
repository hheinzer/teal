#include <math.h>

#include "teal.h"

const double rho = 1.4, vmag = 0.5, alpha = 3 * PI / 180, p = 1;

int main(int argc, char **argv) {
    teal_init(argc, argv);

    const double state[] = {rho, vmag * cos(alpha), vmag * sin(alpha), p};
    Mesh mesh = airfoil_mesh("dat/naca2312.dat", 200, 200, 100, 20);
    mesh_set_boundary_condition(&mesh, "airfoil", "wall", 0, 0);
    mesh_set_boundary_condition(&mesh, "farfield", "characteristic", state, 0);
    mesh_generate(&mesh);
    mesh_print(&mesh);

    Equations eqns = euler_create(&mesh, 0);
    equations_set_initial_state(&eqns, 4, (long[]){D, U, V, P}, state);
    euler_compute_conserved(&eqns);
    euler_print(&eqns);

    Simulation sim = simulation_create(argv[0], &eqns);
    sim.abort_variable = DU;
    sim.abort_residual = 1e-4;
    sim.output_iter = 10000;
    sim.max_iter = sim.output_iter * 10;
    simulation_print(&sim);
    simulation_run(&sim);
    airfoil_polar(&sim, 0, P, state, 1);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);
    teal_finalize();
}
