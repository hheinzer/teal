#include "riemann.h"
#include "teal.h"

static const double dl = 1, ul = 0, pl = 1;
static const double dr = 0.125, ur = 0, pr = 0.1;
static const double x0 = 0.5, max_time = 0.25;
static Function exact;

int main(int argc, char **argv)
{
    teal_init(argc, argv);

    Mesh mesh = mesh_create((double[]){0, 0}, (double[]){1, 1}, (long[]){100, 1});
    mesh_set_periodic_condition(&mesh, "bottom", "top");
    mesh_set_boundary_condition(&mesh, "left", "outflow", 0, 0);
    mesh_set_boundary_condition(&mesh, "right", "outflow", 0, 0);
    mesh_generate(&mesh);
    mesh_print(&mesh);

    char *name[] = {"exact density", "exact velocity-x", "exact velocity-y", "exact pressure"};
    Fields user = {.n_fields = 4, .name = name, .compute = exact};
    Equations eqns = euler_create(&mesh, &user);
    equations_set_initial_condition(&eqns, exact);
    euler_compute_conserved(&eqns);
    euler_print(&eqns);

    Simulation sim = simulation_create(argv[0], &eqns);
    sim.max_time = max_time;
    simulation_print(&sim);
    simulation_run(&sim);
    simulation_error(&sim, D, 0);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);
    teal_finalize();
}

static void exact(const double *x, const double time, const double *, double *u)
{
    riemann(1.4, dl, ul, pl, dr, ur, pr, (x[0] - x0) / time, &u[D], &u[U], &u[P]);
    u[V] = 0;
}
