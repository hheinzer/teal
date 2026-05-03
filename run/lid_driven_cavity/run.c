#include <math.h>

#include "navierstokes.h"

static const double mach = 0.1, reynolds = 100;

static void moving_wall(void *outer_, const double *property, Vector center, double time,
                        const void *context)
{
    NavierStokesPrimitive *outer = outer_;
    const NavierStokesPrimitive *inner = context;
    outer->density = inner->density;
    outer->velocity.x = 2 - inner->velocity.x;
    outer->velocity.y = -inner->velocity.y;
    outer->velocity.z = -inner->velocity.z;
    outer->pressure = inner->pressure;
}

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);

    Vector min_coord = {.x = 0, .y = 0};
    Vector max_coord = {.x = 1, .y = 1};
    Triple num_cells = {.x = 100, .y = 100};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, (Triple){0});

    for (long i = 0; i < mesh->nodes.num; i++) {
        Vector *coord = &mesh->nodes.coord[i];
        double beta = 2;
        coord->x = (1 + (tanh(beta * ((2 * coord->x) - 1)) / tanh(beta))) / 2;
        coord->y = (1 + (tanh(beta * ((2 * coord->y) - 1)) / tanh(beta))) / 2;
    }

    mesh_generate(mesh);
    mesh_summary(mesh);

    NavierStokesPrimitive state = {.density = 1, .pressure = 1 / (1.4 * sq(mach))};

    Equations *eqns = navierstokes_create(mesh);
    equations_set_limiter(eqns, "venkatakrishnan", 1);
    equations_set_boundary_condition(eqns, "left", "wall", 0, 0);
    equations_set_boundary_condition(eqns, "right", "wall", 0, 0);
    equations_set_boundary_condition(eqns, "bottom", "wall", 0, 0);
    equations_set_boundary_condition(eqns, "top", "custom", moving_wall, 0);
    equations_set_initial_condition(eqns, "domain", 0, &state);
    equations_set_property(eqns, "dynamic viscosity", 1 / reynolds);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_iter(sim, 10000);
    simulation_set_out_iter(sim, 1000);
    simulation_set_advance(sim, "implicit euler", 100, 0);
    simulation_summary(sim);

    simulation_run(sim);

    simulation_destroy(sim);
    equations_destroy(eqns);
    mesh_destroy(mesh);

    teal_deinit();
}
