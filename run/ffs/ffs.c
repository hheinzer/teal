#include "teal.h"

int main(int argc, char **argv)
{
    teal_init(argc, argv);

    const double state[] = {1.4, 3, 0, 1};
    Mesh mesh = mesh_read("run/ffs/ffs.geo");
    mesh_set_boundary_condition(&mesh, "inflow", "inflow", state, 0);
    mesh_set_boundary_condition(&mesh, "outflow", "outflow", 0, 0);
    mesh_set_boundary_condition(&mesh, "symmetry", "symmetry", 0, 0);
    mesh_set_boundary_condition(&mesh, "wall", "wall", 0, 0);
    mesh_generate(&mesh);
    mesh_print(&mesh);

    Equations eqns = euler_create(&mesh, 0);
    eqns.flux = euler_flux("hll");
    equations_set_initial_state(&eqns, 4, (long[]){D, U, V, P}, state);
    euler_compute_conserved(&eqns);
    euler_print(&eqns);

    Simulation sim = simulation_create(argv[0], &eqns);
    sim.max_time = 4;
    simulation_print(&sim);
    simulation_run(&sim);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);
    teal_finalize();
}
