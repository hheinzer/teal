#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "euler.h"

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    if (argc > 1 && !strcmp(argv[1], "-h"))
        teal_exit("usage: %s -- [-h] [geo file] [mach number] [alpha]", argv[0]);

    Mesh *mesh = mesh_read((argc > 1 ? argv[1] : "run/airfoil/naca2312.geo"));
    mesh_generate(&mesh);
    mesh_print(mesh);

    const double vmag = (argc > 2 ? strtod(argv[2], 0) : 0.5);
    const double alpha = (argc > 3 ? strtol(argv[3], 0, 0) : 3) * PI / 180;
    const double state[N_VARS] = {
        [D] = 1.4, [U] = vmag * cos(alpha), [V] = vmag * sin(alpha), [P] = 1};
    Equations *eqns = euler_create(mesh);
    equations_set_limiter(eqns, "venk", 1);
    equations_set_boundary_condition(eqns, "airfoil", "slipwall", state, 0);
    equations_set_boundary_condition(eqns, "farfield", "farfield", state, 0);
    equations_set_initial_state(eqns, state);
    equations_print(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_iter(sim, 10000);
    simulation_set_output_iter(sim, 100);
    simulation_set_abort(sim, D, 1e-5);
    simulation_set_advance(sim, "implicit euler");
    simulation_set_cfl(sim, 100);
    simulation_print(sim);
    simulation_run(sim);
    euler_polar(sim, "airfoil", state, 1);

    simulation_free(&sim);
    equations_free(&eqns);
    mesh_free(&mesh);

    teal_finalize();
}
