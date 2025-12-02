#define _GNU_SOURCE
#include <math.h>
#include <stdlib.h>

#include "euler.h"

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    const char *fname = (argc > 1) ? argv[1] : "run/airfoil/naca2312.msh";
    scalar mach = (argc > 2) ? strtod(argv[2], 0) : 0.5;
    scalar alpha = ((argc > 3) ? strtod(argv[3], 0) : 3) * M_PI / 180;

    Mesh *mesh = mesh_read(fname);
    mesh_generate(mesh);
    mesh_summary(mesh);

    Euler farfield = {
        .density = 1.4,
        .velocity = {.x = mach * cos(alpha), .y = mach * sin(alpha)},
        .pressure = 1,
    };

    Equations *eqns = euler_create(mesh);
    equations_set_limiter(eqns, "venkatakrishnan", 1);
    equations_set_boundary_condition(eqns, "airfoil", "slipwall", 0, 0);
    equations_set_boundary_condition(eqns, "farfield", "farfield", &farfield, 0);
    equations_set_initial_state(eqns, "domain", &farfield);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_iter(sim, 10000);
    simulation_set_out_iter(sim, 100);
    simulation_set_termination(sim, "density", 1e-5);
    simulation_set_advance(sim, "implicit euler", 100, 0);
    simulation_summary(sim);

    scalar time = simulation_run(sim);
    euler_polar(sim, "airfoil", &farfield, 1, time);

    teal_finalize();
}
