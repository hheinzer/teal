#include <math.h>

#include "teal.h"

const double alpha = 30 * PI / 180;
const double state[2][4] = {{1.4, 0, 0, 1}, {8, 7.1447096, -4.125, 116.5}};
const double Wx = 8.660254, Wy = -5;
static Function bc;

int main(int argc, char **argv) {
    teal_init(argc, argv);

    Mesh mesh = mesh_create((double[]){0, 0}, (double[]){3.25, 1}, (long[]){520, 160});
    mesh_split_entity(&mesh, "bottom", 0, 1.0 / 6);
    mesh_set_boundary_condition(&mesh, "bottom0", "outflow", 0, 0);
    mesh_set_boundary_condition(&mesh, "bottom1", "wall", 0, 0);
    mesh_set_boundary_condition(&mesh, "right", "outflow", 0, 0);
    mesh_set_boundary_condition(&mesh, "top", "custom", 0, bc);
    mesh_set_boundary_condition(&mesh, "left", "inflow", state[1], 0);
    mesh_generate(&mesh);
    mesh_print(&mesh);

    Equations eqns = euler_create(&mesh, 0);
    equations_set_initial_condition(&eqns, bc);
    euler_compute_conserved(&eqns);
    euler_print(&eqns);

    Simulation sim = simulation_create(argv[0], &eqns);
    sim.max_time = 0.2;
    simulation_print(&sim);
    simulation_run(&sim);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);
    teal_finalize();
}

static void bc(const double *x, const double time, const double *, double *u) {
    const long s = (atan2(x[0] - 1.0 / 6 - Wx * time, x[1] - Wy * time) < alpha);
    u[D] = state[s][0];
    u[U] = state[s][1];
    u[V] = state[s][2];
    u[P] = state[s][3];
}
