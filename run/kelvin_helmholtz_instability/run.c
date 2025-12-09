#define _GNU_SOURCE
#include <math.h>

#include "navier_stokes.h"
#include "teal/utils.h"

scalar viscosity = 0;
Compute initial;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    vector min_coord = {.x = 0, .y = 0};
    vector max_coord = {.x = 1, .y = 1};
    tuple num_cells = {.x = 256, .y = 256};
    bool periodic[] = {true, true};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, periodic);
    mesh_generate(mesh);
    mesh_summary(mesh);

    Equations *eqns = navier_stokes_create(mesh);
    equations_set_limiter(eqns, "venkatakrishnan", 1);
    equations_set_initial_condition(eqns, "domain", initial, 0);
    equations_set_property(eqns, "dynamic viscosity", viscosity);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 1);
    simulation_set_out_time(sim, 0.1);
    simulation_summary(sim);

    simulation_run(sim);

    teal_finalize();
}

void initial(void *variable_, const scalar *property, vector center, scalar time, const void *ctx_)
{
    NavierStokes *variable = variable_;
    if (0.25 < center.y && center.y < 0.75) {
        variable->density = 2;
        variable->velocity.x = 0.5;
    }
    else {
        variable->density = 1;
        variable->velocity.x = -0.5;
    }
    scalar w0 = 0.1, sigma2 = sq(0.05 / 2);
    variable->velocity.y =
        w0 * sin(4 * M_PI * center.x) *
        (exp(-sq(center.y - 0.25) / (2 * sigma2)) + exp(-sq(center.y - 0.75) / (2 * sigma2)));
    variable->pressure = 2.5;
}
