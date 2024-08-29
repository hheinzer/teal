#include <math.h>

#include "euler.h"

Compute rmi;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    Mesh *mesh = mesh_create((Vector3d){0, 0, 0}, (Vector3d){60, 15, 15}, (Vector3l){400, 100, 1});
    mesh_periodic(mesh, "back", "front");
    mesh_generate(&mesh);
    mesh_print(mesh);

    Equations *eqns = euler_create(mesh);
    equations_set_limiter(eqns, "venk", 1);
    equations_set_boundary_condition(eqns, "left", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "right", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "bottom", "symmetry", 0, 0);
    equations_set_boundary_condition(eqns, "top", "symmetry", 0, 0);
    equations_set_initial_condition(eqns, rmi);
    equations_print(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 200);
    simulation_set_output_time(sim, 2);
    simulation_print(sim);
    simulation_run(sim);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);

    teal_finalize();
}

void rmi(double *u, const double *, const double *x, double)
{
    u[P] = 1;
    if (x[X] >= 18 + 2 * cos(2 * PI * x[Y] / 15))
        u[D] = 0.25;
    else
        u[D] = 1;
    if (2 <= x[X] && x[X] <= 6) {
        u[D] = 4.22;
        u[P] = 4.9;
    }
}
