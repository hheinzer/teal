#include "euler.h"

#include <math.h>

#include "mesh.h"
#include "simulation.h"
#include "teal.h"

#define PI 3.14159265358979323846
const double gamma = 1.4, a = 0.1, b = PI, c = 2 * PI;
Function exact, source;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    const double x0[] = {-1, -1, -1}, x1[] = {1, 1, 1};
    const long n_cells[] = {50, 50, 50};
    Mesh mesh = mesh_create(x0, x1, n_cells);
    // Mesh mesh = mesh_read("run/sinewave/tets.geo");
    mesh_periodic(&mesh, "left", "right");
    mesh_periodic(&mesh, "bottom", "top");
    mesh_periodic(&mesh, "back", "front");
    mesh_generate(&mesh);
    mesh_print(&mesh);

    Equations eqns = euler_create(&mesh, 2);
    equations_create_exact(&eqns, exact, 5);
    equations_set_scalar(&eqns, GAMMA, gamma);
    equations_set_source(&eqns, source);
    equations_set_initial_condition(&eqns, exact);
    equations_print(&eqns);

    Simulation sim = simulation_create(&eqns, argv[0]);
    simulation_set_max_time(&sim, 0.5);
    simulation_print(&sim);
    simulation_run(&sim);
    simulation_error(&sim, D, 0);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);

    teal_finalize();
}

void exact(double *u, const double *, const double *x, double time)
{
    u[D] = a * sin(b * (x[X] + x[Y] + x[Z]) - c * time) + 2;
    u[U] = u[V] = u[W] = 1;
    u[P] = (gamma - 1) *
           ((u[D] * u[D]) - 0.5 * u[D] * ((u[U] * u[U]) + (u[V] * u[V]) + (u[W] * u[W])));
}

void source(double *q, const double *, const double *x, double time)
{
    q[D] = a * (3 * b - c) * cos(b * (x[X] + x[Y] + x[Z]) - c * time);
    q[DU] = q[DV] = q[DW] =
        (1.0 / 2.0) * a *
        (b * (gamma - 1) * (4 * a * sin(b * (x[X] + x[Y] + x[Z]) - c * time) + 5) + 6 * b - 2 * c) *
        cos(b * (x[X] + x[Y] + x[Z]) - c * time);
    q[DE] = 0.5 * a *
            (12 * a * b * gamma * sin(b * x[X] + b * x[Y] + b * x[Z] - c * time) -
             4 * a * c * sin(b * x[X] + b * x[Y] + b * x[Z] - c * time) + 15 * b * gamma + 9 * b -
             8 * c) *
            cos(b * x[X] + b * x[Y] + b * x[Z] - c * time);
}
