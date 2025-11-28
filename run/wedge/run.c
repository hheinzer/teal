#define _GNU_SOURCE
#include <math.h>

#include "euler.h"

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    Mesh *mesh = mesh_read("run/wedge/mesh.msh");
    mesh_generate(mesh);
    mesh_summary(mesh);

    scalar alpha = 5 * M_PI / 180;
    scalar velocity = 3;
    Euler inlet = {
        .density = 1.4,
        .velocity = {.x = cos(alpha) * velocity, .y = sin(alpha) * velocity},
        .pressure = 1,
    };

    Equations *eqns = euler_create(mesh);
    equations_set_boundary_condition(eqns, "inlet", "supersonic inflow", &inlet, 0);
    equations_set_boundary_condition(eqns, "outlet", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "airfoil", "wall", 0, 0);
    equations_set_initial_state(eqns, "domain", &inlet);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_iter(sim, 5000);
    simulation_set_out_iter(sim, 50);
    simulation_summary(sim);

    scalar time = simulation_run(sim);
    euler_polar(sim, "airfoil", &inlet, 2, time);

    teal_finalize();
}
