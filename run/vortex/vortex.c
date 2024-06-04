#include <math.h>

#include "teal.h"

const double gamma = 1.4, x0 = 5, y0 = 5, beta = 5, u0 = 1, v0 = 1;
static Function exact;

int main(int argc, char **argv)
{
    teal_init(argc, argv);

    Mesh mesh = mesh_create((double[]){0, 0}, (double[]){10, 10}, (long[]){100, 100});
    mesh_set_periodic_condition(&mesh, "bottom", "top");
    mesh_set_periodic_condition(&mesh, "left", "right");
    mesh_generate(&mesh);
    mesh_print(&mesh);

    char *name[] = {"exact density", "exact velocity-x", "exact velocity-y", "exact pressure"};
    Fields user = {.n_fields = 4, .name = name, .compute = exact};
    Equations eqns = euler_create(&mesh, &user);
    equations_set_initial_condition(&eqns, exact);
    euler_compute_conserved(&eqns);
    euler_print(&eqns);

    Simulation sim = simulation_create(argv[0], &eqns);
    sim.max_time = 10;
    simulation_print(&sim);
    simulation_run(&sim);
    simulation_error(&sim, D, 0);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);
    teal_finalize();
}

static double fwrap(const double x, const double xmin, const double xmax)
{
    if (xmin > xmax) return fwrap(x, xmax, xmin);
    return (x >= 0 ? xmin : xmax) + fmod(x, xmax - xmin);
}

static void exact(const double *x, const double time, const double *, double *u)
{
    const double xi = fwrap(x[X] - u0 * time, 0, 10);
    const double yi = fwrap(x[Y] - v0 * time, 0, 10);
    const double r2 = (xi - x0) * (xi - x0) + (yi - y0) * (yi - y0);
    u[D] =
        pow(1 - (gamma - 1) * beta * beta / (8 * gamma * PI * PI) * exp(1 - r2), 1 / (gamma - 1));
    u[U] = u0 - beta / (2 * PI) * exp(0.5 * (1 - r2)) * (yi - y0);
    u[V] = v0 + beta / (2 * PI) * exp(0.5 * (1 - r2)) * (xi - x0);
    u[P] = pow(u[D], gamma);
}
