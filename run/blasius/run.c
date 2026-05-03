#include <math.h>

#include "navierstokes.h"

static const double mach = 0.1, reynolds = 1000;

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);

    Vector min_coord = {.x = -0.5, .y = 0};
    Vector max_coord = {.x = 1, .y = 0.5};
    Triple num_cells = {.x = 75, .y = 75};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, (Triple){0});

    for (long i = 0; i < mesh->nodes.num; i++) {
        Vector *coord = &mesh->nodes.coord[i];
        double beta = 3;
        coord->x = copysign(expm1(beta * fabs(coord->x)) / expm1(beta), coord->x);
        coord->y = expm1(beta * coord->y) / expm1(beta);
    }

    Vector root = {.x = 0, .y = 0};
    Vector normal = {.x = 1, .y = 0};
    mesh_split(mesh, "bottom", root, normal);

    mesh_generate(mesh);
    mesh_summary(mesh);

    NavierStokesPrimitive farfield = {
        .density = 1, .velocity = {.x = 1}, .pressure = 1 / (1.4 * sq(mach))};

    Equations *eqns = navierstokes_create(mesh);
    equations_set_limiter(eqns, "venkatakrishnan", 1);
    equations_set_boundary_condition(eqns, "left", "subsonic inflow", 0, &farfield);
    equations_set_boundary_condition(eqns, "right", "pressure outflow", 0, &farfield);
    equations_set_boundary_condition(eqns, "bottom-a", "slipwall", 0, 0);
    equations_set_boundary_condition(eqns, "bottom-b", "adiabatic wall", 0, 0);
    equations_set_boundary_condition(eqns, "top", "pressure outflow", 0, &farfield);
    equations_set_initial_condition(eqns, "domain", 0, &farfield);
    equations_set_property(eqns, "dynamic viscosity", 1 / reynolds);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_iter(sim, 10000);
    simulation_set_out_iter(sim, 100);
    simulation_set_termination(sim, "momentum-x", 1e-5);
    simulation_set_advance(sim, "implicit euler", 100, 0);
    simulation_summary(sim);

    simulation_run(sim);

    simulation_destroy(sim);
    equations_destroy(eqns);
    mesh_destroy(mesh);

    teal_deinit();
}
