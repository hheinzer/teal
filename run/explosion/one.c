#include <stdio.h>

#include "euler.h"

const double alpha = 2;
const double ui[N_VARS] = {[D] = 1, [P] = 1};
const double uo[N_VARS] = {[D] = 0.125, [P] = 0.1};
Compute initial, source;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    Mesh *mesh = mesh_create((Vector3d){0, 0, 0}, (Vector3d){1, 1, 1}, (Vector3l){1000, 1, 1});
    mesh_periodic(mesh, "bottom", "top");
    mesh_periodic(mesh, "back", "front");
    mesh_generate(&mesh);
    mesh_print(mesh);

    Equations *eqns = euler_create(mesh);
    equations_set_source(eqns, source, 0);
    equations_set_boundary_condition(eqns, "left", "symmetry", 0, 0);
    equations_set_boundary_condition(eqns, "right", "supersonic outflow", 0, 0);
    equations_set_initial_condition(eqns, initial);
    equations_print(eqns);

    String prefix;
    sprintf(prefix, "%s_alpha_%g", argv[0], alpha);
    Simulation *sim = simulation_create(eqns, prefix);
    simulation_set_max_time(sim, 0.25);
    simulation_print(sim);
    simulation_run(sim);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);

    teal_finalize();
}

void initial(double *u, const double *, const Vector3d x, double)
{
    for (long v = 0; v < N_VARS; ++v) u[v] = (x[X] <= 0.4 ? ui[v] : uo[v]);
}

void source(double *q, const double *u, const Vector3d x, double)
{
    const double fac = -alpha / x[X];
    q[D] = fac * u[DU];
    q[DU] = fac * u[DU] * u[U];
    q[DE] = fac * (u[DE] + u[P]) * u[U];
}
