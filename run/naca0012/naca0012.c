#include <math.h>

#include "airfoil.h"
#include "euler.h"
#include "mesh.h"
#include "simulation.h"
#include "teal.h"

#define PI 3.14159265358979323846
const double rho = 1.4, vmag = 0.8, alpha = 1.25 * PI / 180, p = 1;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    Mesh mesh = mesh_read("run/naca0012/naca0012.geo");
    mesh_generate(&mesh);
    mesh_print(&mesh);

    const double state[N_VARS] = {
        [D] = rho, [U] = vmag * cos(alpha), [V] = vmag * sin(alpha), [P] = p};
    Equations eqns = euler_create(&mesh, 2);
    equations_set_limiter(&eqns, "venk", 1);
    equations_set_initial_state(&eqns, state);
    equations_set_boundary_condition(&eqns, "airfoil", "slipwall", state, 0);
    equations_set_boundary_condition(&eqns, "farfield", "farfield", state, 0);
    equations_print(&eqns);

    Simulation sim = simulation_create(&eqns, argv[0]);
    simulation_set_max_iter(&sim, 100000);
    simulation_set_output_iter(&sim, 10000);
    simulation_set_abort(&sim, D, 1e-5);
    simulation_print(&sim);
    simulation_run(&sim);
    airfoil_polar(&sim, "airfoil", state, 1, D, U, V, W, P);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);

    teal_finalize();
}
