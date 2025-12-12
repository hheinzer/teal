#include <math.h>

#include "equations.h"
#include "navier_stokes.h"
#include "teal/utils.h"

static scalar mach = 0.1, reynolds = 100;
NavierStokes average;
Source source;
Prepare prepare;

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    vector min_coord = {.x = 0, .y = 0};
    vector max_coord = {.x = 9, .y = 3};
    tuple num_cells = {.x = 150, .y = 50};
    bool periodic[] = {true, false};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, periodic);

    for (long i = 0; i < mesh->nodes.num; i++) {
        vector *coord = &mesh->nodes.coord[i];
        scalar a = 4.5, b = 3.5, c = 1.0 / 6;
        coord->y += c * (3 - coord->y) * (1 + tanh(b * (fabs(coord->x - a) - b)));
    }

    mesh_generate(mesh);
    mesh_summary(mesh);

    NavierStokes initial = {.density = 1, .velocity = {.x = 1}, .pressure = 1 / (1.4 * sq(mach))};

    Equations *eqns = navier_stokes_create(mesh);
    equations_set_limiter(eqns, "venkatakrishnan", 1);
    equations_set_boundary_condition(eqns, "bottom", "wall", 0, 0);
    equations_set_boundary_condition(eqns, "top", "wall", 0, 0);
    equations_set_initial_state(eqns, "domain", &initial);
    equations_set_user_source(eqns, source, prepare);
    equations_set_property(eqns, "dynamic viscosity", 1 / reynolds);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 100);
    simulation_set_out_time(sim, 10);
    simulation_set_advance(sim, "implicit euler", 100, 0);
    simulation_summary(sim);

    simulation_run(sim);

    teal_finalize();
}

void prepare(const Equations *eqns, const void *variable_)
{
    equations_average(eqns, "domain", variable_, &average);
}

void source(void *source_, const void *variable_, const scalar *property, vector center,
            scalar time)
{
    NavierStokes *source = source_;
    const NavierStokes *variable = variable_;
    scalar s1 = variable->density * (1 - average.momentum.x / average.density) / 0.001;
    source->momentum.x = s1;
    source->energy = variable->velocity.x * s1;
}
