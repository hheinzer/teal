#include <math.h>

#include "navierstokes.h"

const double Reyn = 250, Vr = 0;
Compute rotating_wall;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    Mesh *mesh = mesh_read("run/cylinder/cylinder.geo");
    mesh_generate(&mesh);
    mesh_print(mesh);

    const double state[N_VARS] = {[D] = 1.4, [U] = 0.1, [P] = 1};
    Equations *eqns = navierstokes_create(mesh);
    equations_set_scalar(eqns, MU, state[D] * state[U] / Reyn);
    equations_set_limiter(eqns, "venk", 1);
    equations_set_boundary_condition(eqns, "wall", "custom", 0, rotating_wall);
    equations_set_boundary_condition(eqns, "farfield", "farfield", state, 0);
    equations_set_initial_state(eqns, state);
    equations_print(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 1000);
    simulation_set_output_time(sim, 10);
    simulation_print(sim);
    simulation_run(sim);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);

    teal_finalize();
}

void rotating_wall(double *ug, const double *ui, const Vector3d x, double)
{
    const double theta = atan2(x[Y], x[X]) - PI / 2;
    const double uw = Vr * cos(theta);
    const double vw = Vr * sin(theta);
    const double ww = 0;
    ug[D] = ui[D];
    ug[U] = 2 * uw - ui[U];
    ug[V] = 2 * vw - ui[V];
    ug[W] = 2 * ww - ui[W];
    ug[P] = ui[P];
}
