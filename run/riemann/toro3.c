#include "euler.h"

Euler left = {.density = 1, .pressure = 1000};
Euler right = {.density = 1, .pressure = 0.01};
scalar x0 = 0.5, max_time = 0.012;
Compute exact;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    vector min_coord = {.x = 0};
    vector max_coord = {.x = 1};
    tuple num_cells = {.x = 100};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, 0);
    mesh_generate(mesh);
    mesh_summary(mesh);

    Equations *eqns = euler_create(mesh);
    equations_create_exact_solution(eqns, exact);
    equations_set_boundary_condition(eqns, "left", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "right", "supersonic outflow", 0, 0);
    equations_set_initial_condition(eqns, "domain", exact, 0);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, max_time);
    simulation_summary(sim);

    scalar time = simulation_run(sim);
    simulation_error(sim, time, 0);

    teal_finalize();
}

void exact(void *variable_, const scalar *property, vector center, scalar time, const void *ctx_)
{
    Euler *variable = variable_;
    scalar gamma = property[EULER_HEAT_CAPACITY_RATIO];
    scalar location = (center.x - x0) / time;
    *variable = euler_riemann(&left, &right, gamma, location);
}
