#include <math.h>

#include "mesh.h"
#include "navierstokes.h"
#include "simulation.h"
#include "teal.h"

const double Reyn = 1000, x = 2;
Modify grading;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    const double x0[] = {-1, 0, 0}, x1[] = {5, 5, 1}, root[] = {0, 0, 0}, normal[] = {-1, 0, 0};
    const long n_cells[] = {50, 100, 1};
    Mesh mesh = mesh_create(x0, x1, n_cells);
    mesh_modify(&mesh, grading);
    mesh_split(&mesh, "bottom", root, normal);
    mesh_periodic(&mesh, "back", "front");
    mesh_generate(&mesh);
    mesh_print(&mesh);

    const double state[N_VARS] = {[D] = 1.4, [U] = 0.1, [P] = 1};
    Equations eqns = navierstokes_create(&mesh, 2);
    equations_set_scalar(&eqns, MU, state[D] * state[U] * x / Reyn);
    equations_set_limiter(&eqns, "venk", 1);
    equations_set_initial_state(&eqns, state);
    equations_set_boundary_condition(&eqns, "left", "subsonic inflow", state, 0);
    equations_set_boundary_condition(&eqns, "right", "subsonic outflow", state, 0);
    equations_set_boundary_condition(&eqns, "bottom-0", "slipwall", 0, 0);
    equations_set_boundary_condition(&eqns, "bottom-1", "wall", 0, 0);
    equations_set_boundary_condition(&eqns, "top", "subsonic outflow", state, 0);
    equations_print(&eqns);

    Simulation sim = simulation_create(&eqns, argv[0]);
    simulation_set_implicit(&sim, 1000, 20, 0.1, 0.9);
    simulation_set_cfl(&sim, 100);
    simulation_set_output_iter(&sim, 10);
    simulation_set_abort(&sim, DU, 1e-5);
    simulation_print(&sim);
    simulation_run(&sim);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);

    teal_finalize();
}

void grading(double *x)
{
    x[Y] = 5 * pow(x[Y] / 5, 1.2);
}
