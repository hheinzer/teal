#include "euler.h"

#include "airfoil.h"
#include "mesh.h"
#include "simulation.h"
#include "teal.h"

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    Mesh mesh = mesh_read("run/cylinder/cylinder.geo");
    mesh_generate(&mesh);
    mesh_print(&mesh);

    const double state[N_VARS] = {[D] = 1.4, [U] = 0.1, [P] = 1};
    Equations eqns = euler_create(&mesh, 2);
    equations_set_limiter(&eqns, "venk", 1);
    equations_set_initial_state(&eqns, state);
    equations_set_boundary_condition(&eqns, "wall", "slipwall", state, 0);
    equations_set_boundary_condition(&eqns, "farfield", "farfield", state, 0);
    equations_print(&eqns);

    Simulation sim = simulation_create(&eqns, argv[0]);
    simulation_set_output_iter(&sim, 1000);
    simulation_set_abort(&sim, D, 3e-5);
    simulation_print(&sim);
    simulation_run(&sim);
    airfoil_polar(&sim, "wall", state, 2, D, U, V, W, P);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);

    teal_finalize();
}
