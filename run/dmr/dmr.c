#include <math.h>

#include "euler.h"

const double alpha = 30 * PI / 180;
const double state0[N_VARS] = {[D] = 1.4, [P] = 1};
const double state1[N_VARS] = {[D] = 8, [U] = 7.1447096, [V] = -4.125, [P] = 116.5};
const double Wx = 8.660254, Wy = -5;
Compute dmr;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    Mesh *mesh = mesh_create((Vector3d){0, 0, 0}, (Vector3d){3.25, 1, 1}, (Vector3l){520, 160, 1});
    mesh_split(mesh, "bottom", (Vector3d){1.0 / 6, 0, 0}, (Vector3d){-1, 0, 0});
    mesh_periodic(mesh, "back", "front");
    mesh_generate(&mesh);
    mesh_print(mesh);

    Equations *eqns = euler_create(mesh);
    equations_set_boundary_condition(eqns, "bottom-a", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "bottom-b", "slipwall", 0, 0);
    equations_set_boundary_condition(eqns, "right", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "top", "custom", 0, dmr);
    equations_set_boundary_condition(eqns, "left", "supersonic inflow", state1, 0);
    equations_set_initial_condition(eqns, dmr);
    equations_print(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 0.2);
    simulation_print(sim);
    simulation_run(sim);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);

    teal_finalize();
}

void dmr(double *u, const double *, const Vector3d x, double time)
{
    if (atan2(x[X] - 1.0 / 6 - Wx * time, x[Y] - Wy * time) >= alpha)
        for (long v = 0; v < N_VARS; ++v) u[v] = state0[v];
    else
        for (long v = 0; v < N_VARS; ++v) u[v] = state1[v];
}
