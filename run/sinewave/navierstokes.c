#include "navierstokes.h"

#include <math.h>

const double gamma = 1.4, mu = 0.1, prandtl = 0.72;
const double a = 0.1, b = PI, c = 2 * PI;
Compute exact, source;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    Mesh *mesh = mesh_create((Vector3d){-1, -1, -1}, (Vector3d){1, 1, 1}, (Vector3l){50, 50, 50});
    mesh_periodic(mesh, "left", "right");
    mesh_periodic(mesh, "bottom", "top");
    mesh_periodic(mesh, "back", "front");
    mesh_generate(&mesh);
    mesh_print(mesh);

    Equations *eqns = navierstokes_create(mesh);
    equations_create_exact(eqns, exact);
    equations_set_scalar(eqns, GAMMA, gamma);
    equations_set_scalar(eqns, MU, mu);
    equations_set_scalar(eqns, PRANDTL, prandtl);
    equations_set_source(eqns, source, 0);
    equations_set_limiter(eqns, "venk", 1);
    equations_set_initial_condition(eqns, exact);
    equations_print(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 0.5);
    simulation_print(sim);
    simulation_run(sim);
    simulation_error(sim, D);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);

    teal_finalize();
}

void exact(double *u, const double *, const Vector3d x, double time)
{
    u[D] = a * sin(b * (x[X] + x[Y] + x[Z]) - c * time) + 2;
    u[U] = u[V] = u[W] = 1;
    u[P] = (gamma - 1) *
           ((u[D] * u[D]) - 0.5 * u[D] * ((u[U] * u[U]) + (u[V] * u[V]) + (u[W] * u[W])));
}

void source(double *q, const double *, const Vector3d x, double time)
{
    q[D] = a * (3 * b - c) * cos(b * (x[X] + x[Y] + x[Z]) - c * time);
    q[DU] = q[DV] = q[DW] =
        0.5 * a *
        (b * (gamma - 1) * (4 * a * sin(b * (x[X] + x[Y] + x[Z]) - c * time) + 5) + 6 * b - 2 * c) *
        cos(b * (x[X] + x[Y] + x[Z]) - c * time);
    q[DE] = 0.5 * a *
            (6 * a * b * gamma * prandtl *
                 sin(2 * b * x[X] + 2 * b * x[Y] + 2 * b * x[Z] - 2 * c * time) -
             2 * a * c * prandtl * sin(2 * b * x[X] + 2 * b * x[Y] + 2 * b * x[Z] - 2 * c * time) +
             6 * pow(b, 2) * gamma * mu * sin(b * x[X] + b * x[Y] + b * x[Z] - c * time) +
             15 * b * gamma * prandtl * cos(b * x[X] + b * x[Y] + b * x[Z] - c * time) +
             9 * b * prandtl * cos(b * x[X] + b * x[Y] + b * x[Z] - c * time) -
             8 * c * prandtl * cos(b * x[X] + b * x[Y] + b * x[Z] - c * time)) /
            prandtl;
}
