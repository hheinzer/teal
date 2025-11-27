#include <stdio.h>
#include <stdlib.h>

#include "euler.h"

scalar alpha = 0;
Euler inner = {.density = 1, .pressure = 1};
Euler outer = {.density = 0.125, .pressure = 0.1};
Compute initial;
Source source;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    if (argc > 1) {
        alpha = strtod(argv[1], 0);
    }

    vector min_coord = {.x = 0};
    vector max_coord = {.x = 1};
    tuple num_cells = {.x = 1000};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, 0);
    mesh_generate(mesh);
    mesh_summary(mesh);

    Equations *eqns = euler_create(mesh);
    equations_set_user_source(eqns, source);
    equations_set_boundary_condition(eqns, "left", "symmetry", 0, 0);
    equations_set_boundary_condition(eqns, "right", "supersonic outflow", 0, 0);
    equations_set_initial_condition(eqns, "domain", initial, 0);
    equations_summary(eqns);

    char prefix[128];
    sprintf(prefix, "%s_alpha_%g", argv[0], alpha);

    Simulation *sim = simulation_create(eqns, prefix);
    simulation_set_max_time(sim, 0.25);
    simulation_summary(sim);

    simulation_run(sim);

    teal_finalize();
}

void initial(void *variable_, const scalar *property, vector center, scalar time, const void *ctx_)
{
    Euler *variable = variable_;
    *variable = (center.x <= 0.4) ? inner : outer;
}

void source(void *source_, const void *variable_, const scalar *property, vector center,
            scalar time)
{
    scalar *source = source_;
    const Euler *variable = variable_;
    scalar fac = -alpha / center.x;
    source[0] = fac * variable->momentum.x;
    source[1] = fac * variable->momentum.x * variable->velocity.x;
    source[2] = 0;
    source[3] = 0;
    source[4] = fac * (variable->energy + variable->pressure) * variable->velocity.x;
}
