#include <stdio.h>

#include "euler.h"
#include "mesh.h"
#include "simulation.h"
#include "teal.h"

const double alpha = 2;
const double ui[N_VARS] = {[D] = 1, [P] = 1};
const double uo[N_VARS] = {[D] = 0.125, [P] = 0.1};
Function initial, source;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    const double x0[] = {0, 0, 0}, x1[] = {1, 1, 1};
    const long n_cells[] = {1000, 1, 1};
    Mesh mesh = mesh_create(x0, x1, n_cells);
    mesh_periodic(&mesh, "bottom", "top");
    mesh_periodic(&mesh, "back", "front");
    mesh_generate(&mesh);
    mesh_print(&mesh);

    Equations eqns = euler_create(&mesh, 2);
    equations_set_source(&eqns, source);
    equations_set_initial_condition(&eqns, initial);
    equations_set_boundary_condition(&eqns, "left", "symmetry", 0, 0);
    equations_set_boundary_condition(&eqns, "right", "supersonic outflow", 0, 0);
    equations_print(&eqns);

    char prefix[128];
    sprintf(prefix, "%s_alpha_%g", argv[0], alpha);
    Simulation sim = simulation_create(&eqns, prefix);
    simulation_set_max_time(&sim, 0.25);
    simulation_print(&sim);
    simulation_run(&sim);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);

    teal_finalize();
}

void initial(double *u, const double *, const double *x, double)
{
    if (x[X] <= 0.4)
        for (long v = 0; v < N_VARS; ++v) u[v] = ui[v];
    else
        for (long v = 0; v < N_VARS; ++v) u[v] = uo[v];
}

void source(double *q, const double *u, const double *x, double)
{
    const double fac = -alpha / x[X];
    q[D] = fac * u[DU];
    q[DU] = fac * u[DU] * u[U];
    q[DE] = fac * (u[DE] + u[P]) * u[U];
}
