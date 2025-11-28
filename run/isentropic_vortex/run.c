#define _GNU_SOURCE
#include <math.h>

#include "euler.h"
#include "teal/utils.h"

vector position = {.x = 5, .y = 5};
vector velocity = {.x = 1, .y = 1};
scalar vortex_strength = 5;
Compute exact;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    vector min_coord = {.x = 0, .y = 0};
    vector max_coord = {.x = 10, .y = 10};
    tuple num_cells = {.x = 100, .y = 100};
    bool periodic[] = {true, true};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, periodic);
    mesh_generate(mesh);
    mesh_summary(mesh);

    Equations *eqns = euler_create(mesh);
    equations_create_exact_solution(eqns, exact);
    equations_set_initial_condition(eqns, "domain", exact, 0);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 10);
    simulation_set_out_time(sim, 0.1);
    simulation_summary(sim);

    scalar time = simulation_run(sim);
    simulation_error(sim, time, 0);

    teal_finalize();
}

static scalar wrap(scalar x, scalar xmin, scalar xmax)
{
    return (x >= 0 ? xmin : xmax) + fmod(x, xmax - xmin);
}

void exact(void *variable_, const scalar *property, vector center, scalar time, const void *ctx_)
{
    Euler *variable = variable_;
    scalar gamma = property[0];

    scalar dx = wrap(center.x - (velocity.x * time), 0, 10) - position.x;
    scalar dy = wrap(center.y - (velocity.y * time), 0, 10) - position.y;
    scalar r2 = (dx * dx) + (dy * dy);

    variable->density =
        pow(1 - ((gamma - 1) * sq(vortex_strength) / (8 * gamma * sq(M_PI)) * exp(1 - r2)),
            1 / (gamma - 1));
    variable->velocity.x = velocity.x - (vortex_strength / (2 * M_PI) * exp((1 - r2) / 2) * dy);
    variable->velocity.y = velocity.y + (vortex_strength / (2 * M_PI) * exp((1 - r2) / 2) * dx);
    variable->velocity.z = 0;
    variable->pressure = pow(variable->density, gamma);
}
