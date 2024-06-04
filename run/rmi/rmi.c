#include <math.h>

#include "teal.h"

static Function initial;

int main(int argc, char **argv)
{
    teal_init(argc, argv);

    Mesh mesh = mesh_create((double[]){0, 0}, (double[]){60, 15}, (long[]){300, 75});
    mesh_set_boundary_condition(&mesh, "bottom", "symmetry", 0, 0);
    mesh_set_boundary_condition(&mesh, "top", "symmetry", 0, 0);
    mesh_set_boundary_condition(&mesh, "left", "outflow", 0, 0);
    mesh_set_boundary_condition(&mesh, "right", "outflow", 0, 0);
    mesh_generate(&mesh);
    mesh_print(&mesh);

    Equations eqns = euler_create(&mesh, 0);
    equations_set_space_order(&eqns, 2, "venk", 1);
    equations_set_initial_condition(&eqns, initial);
    euler_compute_conserved(&eqns);
    euler_print(&eqns);

    Simulation sim = simulation_create(argv[0], &eqns);
    sim.max_time = 200;
    simulation_print(&sim);
    simulation_run(&sim);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);
    teal_finalize();
}

static void initial(const double *x, const double, const double *, double *u)
{
    u[P] = 1;
    if (x[X] >= 18 + 2 * cos(2 * PI * x[Y] / 15)) {
        u[D] = 0.25;
    }
    else {
        u[D] = 1;
    }
    if (2 <= x[X] && x[X] <= 6) {
        u[D] = 4.22;
        u[P] = 4.9;
    }
}
