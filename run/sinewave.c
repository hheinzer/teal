#include <math.h>

#include "teal.h"

const double gamma = 1.4, a = 0.1, b = PI, c = 2 * PI;
static Function exact, source;

int main(int argc, char **argv) {
    teal_init(argc, argv);

    Mesh mesh = mesh_create((double[]){-1, -1}, (double[]){1, 1}, (long[]){100, 100});
    mesh_set_periodic_condition(&mesh, "left", "right");
    mesh_set_periodic_condition(&mesh, "bottom", "top");
    mesh_generate(&mesh);
    mesh_print(&mesh);

    char *name[] = {"exact density", "exact velocity-x", "exact velocity-y", "exact pressure"};
    Fields user = {.nu = 4, .name = name, .compute = exact};
    Equations eqns = euler_create(&mesh, &user);
    eqns.source = source;
    equations_set_initial_condition(&eqns, exact);
    euler_compute_conserved(&eqns);
    euler_print(&eqns);

    Simulation sim = simulation_create(argv[0], &eqns);
    sim.max_time = c;
    simulation_print(&sim);
    simulation_run(&sim);
    simulation_error(&sim, D, 0);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);
    teal_finalize();
}

static void exact(const double *x, const double time, const double *, double *u) {
    u[D] = 2 + a * sin(b * (x[0] + x[1]) - c * time);
    u[V] = u[U] = 1;
    u[P] = (gamma - 1) * (u[D] * u[D] - 0.5 * u[D] * (u[U] * u[U] + u[V] * u[V]));
}

static void source(const double *x, const double time, const double *, double *q) {
    q[D] = a * (2 * b - c) * cos(b * (x[0] + x[1]) - c * time);
    q[DV] = q[DU] =
        a * (b * (gamma - 1) * (2 * a * sin(b * (x[0] + x[1]) - c * time) + 3) + 2 * b - c) *
        cos(b * (x[0] + x[1]) - c * time);
    q[DE] = 2 * a *
            (2 * a * b * gamma * sin(b * x[0] + b * x[1] - c * time) -
             a * c * sin(b * x[0] + b * x[1] - c * time) + 3 * b * gamma + b - 2 * c) *
            cos(b * x[0] + b * x[1] - c * time);
}
