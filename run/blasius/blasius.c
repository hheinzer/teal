#include <math.h>

#include "navierstokes.h"

const double Reyn = 1000, x0 = 2;
Modify grading;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    Mesh *mesh = mesh_create((Vector3d){-1, 0, 0}, (Vector3d){5, 5, 1}, (Vector3l){150, 150, 1});
    mesh_modify(mesh, grading);
    mesh_split(mesh, "bottom", (Vector3d){0, 0, 0}, (Vector3d){-1, 0, 0});
    mesh_periodic(mesh, "back", "front");
    mesh_generate(&mesh);
    mesh_print(mesh);

    const double state[N_VARS] = {[D] = 1.4, [U] = 0.2, [P] = 1};
    Equations *eqns = navierstokes_create(mesh);
    equations_set_scalar(eqns, MU, state[D] * state[U] * x0 / Reyn);
    equations_set_limiter(eqns, "venk", 1);
    equations_set_initial_state(eqns, state);
    equations_set_boundary_condition(eqns, "left", "farfield", state, 0);
    equations_set_boundary_condition(eqns, "right", "subsonic outflow", state, 0);
    equations_set_boundary_condition(eqns, "bottom-a", "slipwall", 0, 0);
    equations_set_boundary_condition(eqns, "bottom-b", "wall", 0, 0);
    equations_set_boundary_condition(eqns, "top", "farfield", state, 0);
    equations_print(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_output_iter(sim, 50);
    simulation_set_abort(sim, DU, 1e-5);
    simulation_set_advance(sim, "implicit euler");
    simulation_set_cfl(sim, 100);
    simulation_print(sim);
    simulation_run(sim);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);

    teal_finalize();
}

void grading(Vector3d x)
{
    const double r = 3;
    x[Y] = 5 * expm1(r * x[Y] / 5) / expm1(r);
}
