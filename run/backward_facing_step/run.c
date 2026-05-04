#include "navierstokes.h"

static const double mach = 0.1, reynolds = 400;

static void inlet(void *primitive_, const double *property, Vector center, double time,
                  const void *context)
{
    NavierStokesPrimitive *primitive = primitive_;
    double gamma = property[NAVIERSTOKES_HEAT_CAPACITY_RATIO];
    primitive->density = 1;
    primitive->velocity.x = 6 * (center.y - 1) * (2 - center.y);
    primitive->pressure = 1 / (gamma * sq(mach));
}

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);

    Vector min_coord = {.x = 0, .y = 0};
    Vector max_coord = {.x = 30, .y = 2};
    Triple num_cells = {.x = 300, .y = 40};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, (Triple){0});

    Vector root = {.y = 1};
    Vector normal = {.y = 1};
    mesh_split(mesh, "left", root, normal);

    mesh_generate(mesh);
    mesh_summary(mesh);

    NavierStokesPrimitive farfield = {
        .density = 1, .velocity = {.x = 1}, .pressure = 1 / (1.4 * sq(mach))};

    Equations *eqns = navierstokes_create(mesh);
    equations_set_limiter(eqns, "venkatakrishnan", 1);
    equations_set_boundary_condition(eqns, "left-a", "wall", 0, 0);
    equations_set_boundary_condition(eqns, "left-b", "custom", inlet, 0);
    equations_set_boundary_condition(eqns, "right", "pressure outflow", 0, &farfield);
    equations_set_boundary_condition(eqns, "bottom", "wall", 0, 0);
    equations_set_boundary_condition(eqns, "top", "wall", 0, 0);
    equations_set_initial_condition(eqns, "domain", 0, &farfield);
    equations_set_property(eqns, "dynamic viscosity", 2 / reynolds);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_iter(sim, 10000);
    simulation_set_out_iter(sim, 1000);
    simulation_set_termination(sim, "momentum-x", 1e-5);
    simulation_set_advance(sim, "implicit euler", 10, 0);
    simulation_summary(sim);

    simulation_run(sim);

    simulation_destroy(sim);
    equations_destroy(eqns);
    mesh_destroy(mesh);

    teal_deinit();
}
