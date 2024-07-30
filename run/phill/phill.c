#include "mesh.h"
#include "navierstokes.h"
#include "simulation.h"
#include "teal.h"

const double Reyn = 0.72 * 10565;
Prepare prepare;
Function source;
double force[N_DIMS], avg[N_VARS];

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    Mesh mesh = mesh_read("run/phill/phill.geo");
    mesh_periodic(&mesh, "inlet", "outlet");
    mesh_periodic(&mesh, "back", "front");
    mesh_generate(&mesh);
    mesh_print(&mesh);

    const double state[N_VARS] = {[D] = 1.4, [U] = 0.2, [P] = 1};
    Equations eqns = navierstokes_create(&mesh, 2);
    equations_set_scalar(&eqns, MU, state[D] * state[U] / Reyn);
    equations_set_prepare(&eqns, prepare);
    equations_set_source(&eqns, source);
    equations_set_limiter(&eqns, "venk", 1);
    equations_set_initial_state(&eqns, state);
    equations_set_boundary_condition(&eqns, "wall", "wall", 0, 0);
    equations_print(&eqns);

    Simulation sim = simulation_create(&eqns, argv[0]);
    simulation_set_max_time(&sim, 100);
    simulation_set_output_time(&sim, 1);
    simulation_print(&sim);
    simulation_run(&sim);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);

    teal_finalize();
}

void prepare(const Equations *eqns)
{
    navierstokes_body_force(eqns, force);
    equations_average(eqns, avg);
}

void source(double *q, const double *u, const double *, double)
{
    const double tau = 1.0e-3;
    const double S1 = force[X] + u[D] * (0.2 - avg[DU] / avg[D]) / tau;
    q[DU] = S1;
    q[DE] = u[U] * S1;
}
