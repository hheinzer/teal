#include <math.h>

#include "euler.h"

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);

    Mesh *mesh = mesh_read("run/supersonic_wedge/mesh.msh");
    mesh_generate(mesh);
    mesh_summary(mesh);

    double alpha = 5 * M_PI / 180;
    double velocity = 3;
    EulerPrimitive inlet = {
        .density = 1.4,
        .velocity = {.x = cos(alpha) * velocity, .y = sin(alpha) * velocity},
        .pressure = 1,
    };

    Equations *eqns = euler_create(mesh);
    equations_set_boundary_condition(eqns, "inlet", "supersonic inflow", 0, &inlet);
    equations_set_boundary_condition(eqns, "outlet", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "airfoil", "slipwall", 0, 0);
    equations_set_initial_condition(eqns, "domain", 0, &inlet);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_iter(sim, 5000);
    simulation_set_out_iter(sim, 500);
    simulation_summary(sim);

    double time = simulation_run(sim);
    euler_polar(sim, "airfoil", inlet, 2, time);

    simulation_destroy(sim);
    equations_destroy(eqns);
    mesh_destroy(mesh);

    teal_deinit();
}
