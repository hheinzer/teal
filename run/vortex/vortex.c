#include <math.h>

#include "euler.h"
#include "mesh.h"
#include "simulation.h"
#include "teal.h"

#define PI 3.14159265358979323846
const double gamma = 1.4, xx = 5, yy = 5, beta = 5, uu = 1, vv = 1;
Function exact;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    const double x0[] = {0, 0, 0}, x1[] = {10, 10, 10};
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
    equations_set_initial_condition(&eqns, exact);
    equations_print(&eqns);

    Simulation sim = simulation_create(&eqns, argv[0]);
    simulation_set_max_time(&sim, 10);
    simulation_print(&sim);
    simulation_run(&sim);
    simulation_error(&sim, D, 0);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);

    teal_finalize();
}

double fwrap(const double x, const double xmin, const double xmax)
{
    if (xmin > xmax) return fwrap(x, xmax, xmin);
    return (x >= 0 ? xmin : xmax) + fmod(x, xmax - xmin);
}

void exact(double *u, const double *, const double *x, double time)
{
    const double xi = fwrap(x[X] - uu * time, 0, 10);
    const double yi = fwrap(x[Y] - vv * time, 0, 10);
    const double r2 = (xi - xx) * (xi - xx) + (yi - yy) * (yi - yy);
    u[D] =
        pow(1 - (gamma - 1) * beta * beta / (8 * gamma * PI * PI) * exp(1 - r2), 1 / (gamma - 1));
    u[U] = uu - beta / (2 * PI) * exp(0.5 * (1 - r2)) * (yi - yy);
    u[V] = vv + beta / (2 * PI) * exp(0.5 * (1 - r2)) * (xi - xx);
    u[W] = 0;
    u[P] = pow(u[D], gamma);
}
