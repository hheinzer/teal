#include <math.h>

#include "euler.h"

const double ui[N_VARS] = {[D] = 1, [P] = 1};
const double uo[N_VARS] = {[D] = 0.125, [P] = 0.1};
Compute initial;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    Mesh *mesh =
        mesh_create((Vector3d){-1, -1, -1}, (Vector3d){1, 1, 1}, (Vector3l){100, 100, 100});
    mesh_generate(&mesh);
    mesh_print(mesh);

    Equations *eqns = euler_create(mesh);
    equations_set_boundary_condition(eqns, "left", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "right", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "bottom", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "top", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "back", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "front", "supersonic outflow", 0, 0);
    equations_set_initial_condition(eqns, initial);
    equations_print(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 0.25);
    simulation_print(sim);
    simulation_run(sim);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);

    teal_finalize();
}

void initial(double *u, const double *, const double *x, double)
{
    const double r = sqrt((x[X] * x[X]) + (x[Y] * x[Y]) + (x[Z] * x[Z]));
    for (long v = 0; v < N_VARS; ++v) u[v] = (r <= 0.4 ? ui[v] : uo[v]);
}
