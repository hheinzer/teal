#include <stdio.h>
#include <stdlib.h>

#include "euler.h"

static double alpha;

static void initial(void *primitive_, const double *property, Vector center, double time,
                    const void *context)
{
    EulerPrimitive *primitive = primitive_;
    EulerPrimitive inner = {.density = 1, .pressure = 1};
    EulerPrimitive outer = {.density = 0.125, .pressure = 0.1};
    *primitive = (center.x <= 0.4) ? inner : outer;
}

static void source(void *source_, const double *property, Vector center, double time,
                   const void *context)
{
    EulerConserved *source = source_;
    const EulerPrimitive *primitive = context;
    double gamma = property[EULER_HEAT_CAPACITY_RATIO];
    double factor = -alpha / center.x;
    double energy = (primitive->pressure / (gamma - 1)) +
                    (0.5 * primitive->density * sq(primitive->velocity.x));
    source->density = factor * primitive->density * primitive->velocity.x;
    source->momentum.x = factor * primitive->density * sq(primitive->velocity.x);
    source->momentum.y = source->momentum.z = 0;
    source->energy = factor * (energy + primitive->pressure) * primitive->velocity.x;
}

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);

    alpha = (argc > 1) ? strtod(argv[1], 0) : 0;

    Vector min_coord = {0, 0, 0};
    Vector max_coord = {1, 0, 0};
    Triple num_cells = {1000, 0, 0};
    Triple periodic = {0, 0, 0};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, periodic);
    mesh_generate(mesh);
    mesh_summary(mesh);

    Equations *eqns = euler_create(mesh);
    equations_set_boundary_condition(eqns, "left", "symmetry", 0, 0);
    equations_set_boundary_condition(eqns, "right", "supersonic outflow", 0, 0);
    equations_set_initial_condition(eqns, "domain", initial, 0);
    equations_set_source(eqns, source, 0);
    equations_summary(eqns);

    String prefix;
    sprintf(prefix, "%s_alpha_%g", argv[0], alpha);

    Simulation *sim = simulation_create(eqns, prefix);
    simulation_set_max_time(sim, 0.25);
    simulation_summary(sim);

    simulation_run(sim);

    simulation_destroy(sim);
    equations_destroy(eqns);
    mesh_destroy(mesh);

    teal_deinit();
}
