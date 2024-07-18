#include "euler.h"
#include "euler/riemann.h"
#include "mesh.h"
#include "simulation.h"
#include "teal.h"

const double gamma = 1.4;
const double dl = 1, ul = 0, pl = 1;
const double dr = 0.125, ur = 0, pr = 0.1;
const double xd = 0.5, max_time = 0.25;
Function exact;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    const double x0[] = {0, 0, 0}, x1[] = {1, 1, 1};
    const long n_cells[] = {100, 1, 1};
    Mesh mesh = mesh_create(x0, x1, n_cells);
    mesh_periodic(&mesh, "bottom", "top");
    mesh_periodic(&mesh, "back", "front");
    mesh_generate(&mesh);
    mesh_print(&mesh);

    Equations eqns = euler_create(&mesh, 2);
    equations_create_exact(&eqns, exact, 5);
    equations_set_scalar(&eqns, GAMMA, gamma);
    equations_set_initial_condition(&eqns, exact);
    equations_set_boundary_condition(&eqns, "left", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(&eqns, "right", "supersonic outflow", 0, 0);
    equations_print(&eqns);

    Simulation sim = simulation_create(&eqns, argv[0]);
    simulation_set_max_time(&sim, max_time);
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
    riemann(&u[D], &u[U], &u[P], dl, ul, pl, dr, ur, pr, (x[X] - xd) / time, gamma);
    u[V] = u[W] = 0;
}
