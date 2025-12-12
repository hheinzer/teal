#define _GNU_SOURCE
#include "euler.h"

#include <math.h>

#include "teal/utils.h"
#include "teal/vector.h"

scalar a = 0.1, b = M_PI, c = 2 * M_PI;
Compute exact;
Source source;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    vector min_coord = {-1, -1, -1};
    vector max_coord = {1, 1, 1};
    tuple num_cells = {50, 50, 50};
    bool periodic[] = {true, true, true};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, periodic);
    mesh_generate(mesh);
    mesh_summary(mesh);

    Equations *eqns = euler_create(mesh);
    equations_create_exact_solution(eqns, exact);
    equations_set_limiter(eqns, "venkatakrishnan", 1);
    equations_set_initial_condition(eqns, "domain", exact, 0);
    equations_set_user_source(eqns, source, 0);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 0.5);
    simulation_set_out_time(sim, 0.05);
    simulation_summary(sim);

    scalar time = simulation_run(sim);
    simulation_error(sim, time, 0);

    teal_finalize();
}

void exact(void *variable_, const scalar *property, vector center, scalar time, const void *ctx_)
{
    Euler *variable = variable_;
    scalar gamma = property[EULER_HEAT_CAPACITY_RATIO];
    scalar phase = (b * (center.x + center.y + center.z)) - (c * time);
    variable->density = 2 + (a * sin(phase));
    variable->velocity.x = 1;
    variable->velocity.y = 1;
    variable->velocity.z = 1;
    variable->pressure = (gamma - 1) * (sq(variable->density) -
                                        (variable->density * vector_norm2(variable->velocity) / 2));
}

void source(void *source_, const void *variable_, const scalar *property, vector center,
            scalar time)
{
    Euler *source = source_;
    scalar gamma = property[EULER_HEAT_CAPACITY_RATIO];
    scalar phase = (b * (center.x + center.y + center.z)) - (c * time);
    scalar sin_phase = sin(phase);
    scalar cos_phase = cos(phase);
    source->density = a * ((3 * b) - c) * cos_phase;
    source->momentum.x = source->momentum.y = source->momentum.z =
        a * (b * ((gamma - 1) * ((4 * a * sin_phase) + 5)) + (6 * b) - (2 * c)) * cos_phase / 2;
    source->energy = a *
                     (((12 * a * b * gamma * sin_phase) - (4 * a * c * sin_phase)) +
                      (15 * b * gamma) + (9 * b) - (8 * c)) *
                     cos_phase / 2;
}
