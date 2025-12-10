#define _GNU_SOURCE
#include <math.h>

#include "navier_stokes.h"

scalar mach = 0.1, reynolds = 1600;
Compute initial;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    vector min_coord = {-M_PI, -M_PI, -M_PI};
    vector max_coord = {M_PI, M_PI, M_PI};
    tuple num_cells = {64, 64, 64};
    bool periodic[] = {true, true, true};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, periodic);
    mesh_generate(mesh);
    mesh_summary(mesh);

    Equations *eqns = navier_stokes_create(mesh);
    equations_set_limiter(eqns, "venkatakrishnan", 1);
    equations_set_initial_condition(eqns, "domain", initial, 0);
    equations_set_property(eqns, "dynamic viscosity", 1 / reynolds);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 10);
    simulation_set_out_time(sim, 0.1);
    simulation_summary(sim);

    simulation_run(sim);

    teal_finalize();
}

void initial(void *variable_, const scalar *property, vector center, scalar time, const void *ctx_)
{
    NavierStokes *variable = variable_;
    scalar gamma = property[NAVIER_STOKES_HEAT_CAPACITY_RATIO];
    scalar pressure0 = 1 / (gamma * mach * mach);
    variable->velocity.x = sin(center.x) * cos(center.y) * cos(center.z);
    variable->velocity.y = -cos(center.x) * sin(center.y) * cos(center.z);
    variable->velocity.z = 0;
    variable->pressure =
        pressure0 + ((cos(2 * center.x) + cos(2 * center.y)) * (cos(2 * center.z) + 2) / 16);
    variable->density = variable->pressure / pressure0;
}
