#include <math.h>

#include "euler.h"
#include "mesh.h"
#include "simulation.h"
#include "teal.h"

#define PI 3.14159265358979323846
Function initial;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    const double x0[] = {0, 0, 0}, x1[] = {60, 15, 15};
    const long n_cells[] = {400, 100, 1};
    Mesh mesh = mesh_create(x0, x1, n_cells);
    mesh_periodic(&mesh, "back", "front");
    mesh_generate(&mesh);
    mesh_print(&mesh);

    Equations eqns = euler_create(&mesh, 2);
    equations_set_limiter(&eqns, "venk", 1);
    equations_set_initial_condition(&eqns, initial);
    equations_set_boundary_condition(&eqns, "left", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(&eqns, "right", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(&eqns, "bottom", "symmetry", 0, 0);
    equations_set_boundary_condition(&eqns, "top", "symmetry", 0, 0);
    equations_print(&eqns);

    Simulation sim = simulation_create(&eqns, argv[0]);
    simulation_set_max_time(&sim, 200);
    simulation_set_output_iter(&sim, 100);
    simulation_print(&sim);
    simulation_run(&sim);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);
    teal_finalize();
}

void initial(double *u, const double *, const double *x, double)
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
