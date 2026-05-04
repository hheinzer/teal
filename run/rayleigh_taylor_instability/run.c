#include <math.h>

#include "euler.h"

static const double gravity = 0.1, density_high = 2, density_low = 1, pressure_ref = 2.5;

static void initial(void *primitive_, const double *property, Vector center, double time,
                    const void *context)
{
    EulerPrimitive *primitive = primitive_;
    primitive->density = (center.y > 0) ? density_high : density_low;
    primitive->velocity.y =
        0.01 * (1 + cos(4 * M_PI * center.x)) * (1 + cos(3 * M_PI * center.y)) / 4;
    primitive->pressure = pressure_ref - (gravity * primitive->density * center.y);
}

static void source(void *source_, const double *property, Vector center, double time,
                   const void *context)
{
    EulerConserved *source = source_;
    const EulerPrimitive *primitive = context;
    source->momentum.y = -gravity * primitive->density;
    source->energy = -gravity * primitive->density * primitive->velocity.y;
}

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);

    Vector min_coord = {.x = -0.25, .y = -0.75};
    Vector max_coord = {.x = 0.25, .y = 0.75};
    Triple num_cells = {.x = 100, .y = 300};
    Triple periodic = {.x = 1};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, periodic);
    mesh_generate(mesh);
    mesh_summary(mesh);

    Equations *eqns = euler_create(mesh);
    equations_set_limiter(eqns, "venkatakrishnan", 1);
    equations_set_boundary_condition(eqns, "bottom", "symmetry", 0, 0);
    equations_set_boundary_condition(eqns, "top", "symmetry", 0, 0);
    equations_set_initial_condition(eqns, "domain", initial, 0);
    equations_set_source(eqns, source, 0);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 10);
    simulation_set_out_time(sim, 1);
    simulation_summary(sim);

    simulation_run(sim);

    simulation_destroy(sim);
    equations_destroy(eqns);
    mesh_destroy(mesh);

    teal_deinit();
}
