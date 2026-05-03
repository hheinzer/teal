#include <math.h>

#include "navierstokes.h"

static const double mach = 0.1, reynolds = 100;
static NavierStokesPrimitive average;

static void prepare(const Equations *eqns, const void *primitive)
{
    equations_average(eqns, "domain", primitive, &average);
}

static void source(void *source_, const double *property, Vector center, double time,
                   const void *context)
{
    NavierStokesConserved *source = source_;
    const NavierStokesPrimitive *variable = context;
    double s1 = variable->density * (1 - average.velocity.x) / 0.001;
    source->momentum.x = s1;
    source->energy = variable->velocity.x * s1;
}

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);

    Vector min_coord = {.x = 0, .y = 0};
    Vector max_coord = {.x = 9, .y = 3};
    Triple num_cells = {.x = 150, .y = 50};
    Triple periodic = {.x = 1, .y = 0};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, periodic);

    for (long i = 0; i < mesh->nodes.num; i++) {
        Vector *coord = &mesh->nodes.coord[i];
        double a = 4.5, b = 3.5, c = 1.0 / 6;
        coord->y += c * (3 - coord->y) * (1 + tanh(b * (fabs(coord->x - a) - b)));
    }

    mesh_generate(mesh);
    mesh_summary(mesh);

    NavierStokesPrimitive initial = {
        .density = 1, .velocity = {.x = 1}, .pressure = 1 / (1.4 * sq(mach))};

    Equations *eqns = navierstokes_create(mesh);
    equations_set_limiter(eqns, "venkatakrishnan", 1);
    equations_set_boundary_condition(eqns, "bottom", "wall", 0, 0);
    equations_set_boundary_condition(eqns, "top", "wall", 0, 0);
    equations_set_initial_condition(eqns, "domain", 0, &initial);
    equations_set_source(eqns, source, prepare);
    equations_set_property(eqns, "dynamic viscosity", 1 / reynolds);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 100);
    simulation_set_out_time(sim, 10);
    simulation_set_advance(sim, "implicit euler", 100, 0);
    simulation_summary(sim);

    simulation_run(sim);

    simulation_destroy(sim);
    equations_destroy(eqns);
    mesh_destroy(mesh);

    teal_deinit();
}
