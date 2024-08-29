#include "navierstokes.h"

const double Reyn = 100, u0 = 0.2;
Prepare prepare;
Compute source;
double avg[N_VARS];

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    Mesh *mesh = mesh_read("run/phill/phill.geo");
    mesh_periodic(mesh, "inlet", "outlet");
    mesh_periodic(mesh, "back", "front");
    mesh_generate(&mesh);
    mesh_print(mesh);

    const double state[N_VARS] = {[D] = 1.4, [U] = u0, [P] = 1};
    Equations *eqns = navierstokes_create(mesh);
    equations_set_scalar(eqns, MU, state[D] * state[U] / Reyn);
    equations_set_source(eqns, source, prepare);
    equations_set_limiter(eqns, "venk", 1);
    equations_set_boundary_condition(eqns, "wall", "wall", 0, 0);
    equations_set_initial_state(eqns, state);
    equations_print(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 900);
    simulation_set_output_time(sim, 9);
    simulation_set_advance(sim, "implicit euler");
    simulation_set_cfl(sim, 100);
    simulation_print(sim);
    simulation_run(sim);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);

    teal_finalize();
}

void prepare(const Equations *eqns)
{
    equations_average(eqns, avg);
}

void source(double *q, const double *u, const Vector3d, double)
{
    const double tau = 0.01;
    const double S1 = u[D] * (u0 - avg[DU] / avg[D]) / tau;
    q[DU] = S1;
    q[DE] = u[U] * S1;
}
