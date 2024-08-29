#include <math.h>

#include "euler.h"

const double gamma = 1.4;
const double x0 = 5, y0 = 5, u0 = 1, v0 = 1, beta = 5;
Compute exact;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    Mesh *mesh = mesh_create((Vector3d){0, 0, 0}, (Vector3d){10, 10, 10}, (Vector3l){50, 50, 50});
    mesh_periodic(mesh, "left", "right");
    mesh_periodic(mesh, "bottom", "top");
    mesh_periodic(mesh, "back", "front");
    mesh_generate(&mesh);
    mesh_print(mesh);

    Equations *eqns = euler_create(mesh);
    equations_create_exact(eqns, exact);
    equations_set_scalar(eqns, GAMMA, gamma);
    equations_set_initial_condition(eqns, exact);
    equations_print(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 10);
    simulation_print(sim);
    simulation_run(sim);
    simulation_error(sim, D);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);

    teal_finalize();
}

double fwrap(const double x, const double xmin, const double xmax)
{
    return (x >= 0 ? xmin : xmax) + fmod(x, xmax - xmin);
}

void exact(double *u, const double *, const Vector3d x, double time)
{
    const double dx = fwrap(x[X] - u0 * time, 0, 10) - x0;
    const double dy = fwrap(x[Y] - v0 * time, 0, 10) - y0;
    const double r2 = pow(dx, 2) + pow(dy, 2);
    u[D] = pow(1 - (gamma - 1) * pow(beta, 2) / (8 * gamma * pow(PI, 2)) * exp(1 - r2),
               1 / (gamma - 1));
    u[U] = u0 - beta / (2 * PI) * exp((1 - r2) / 2) * dy;
    u[V] = v0 + beta / (2 * PI) * exp((1 - r2) / 2) * dx;
    u[W] = 0;
    u[P] = pow(u[D], gamma);
}
