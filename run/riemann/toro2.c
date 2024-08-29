#include "euler.h"

const double gamma = 1.4;
const double dl = 1, ul = -2, pl = 0.4;
const double dr = 1, ur = 2, pr = 0.4;
const double x0 = 0.5, max_time = 0.15;
Compute exact;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    Mesh *mesh = mesh_create((Vector3d){0, 0, 0}, (Vector3d){1, 1, 1}, (Vector3l){100, 1, 1});
    mesh_periodic(mesh, "bottom", "top");
    mesh_periodic(mesh, "back", "front");
    mesh_generate(&mesh);
    mesh_print(mesh);

    Equations *eqns = euler_create(mesh);
    equations_create_exact(eqns, exact);
    equations_set_scalar(eqns, GAMMA, gamma);
    equations_set_boundary_condition(eqns, "left", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "right", "supersonic outflow", 0, 0);
    equations_set_initial_condition(eqns, exact);
    equations_print(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, max_time);
    simulation_print(sim);
    simulation_run(sim);
    simulation_error(sim, D);

    mesh_write(mesh, argv[0]);
    equations_write(eqns, argv[0], 0, 0);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);

    teal_finalize();
}

void exact(double *u, const double *, const Vector3d x, double time)
{
    euler_riemann(&u[D], &u[U], &u[P], dl, ul, pl, dr, ur, pr, (x[X] - x0) / time, gamma);
    u[V] = u[W] = 0;
}
