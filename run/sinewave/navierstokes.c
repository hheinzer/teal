#include "navierstokes.h"

#include <math.h>

static const double a = 0.1, b = M_PI, c = 2 * M_PI;

static void exact(void *primitive_, const double *property, Vector center, double time,
                  const void *context)
{
    (void)context;
    NavierStokesPrimitive *primitive = primitive_;
    double gamma = property[NAVIERSTOKES_HEAT_CAPACITY_RATIO];
    double phase = (b * (center.x + center.y + center.z)) - (c * time);
    primitive->density = 2 + (a * sin(phase));
    primitive->velocity.x = primitive->velocity.y = primitive->velocity.z = 1;
    primitive->pressure =
        (gamma - 1) *
        (sq(primitive->density) - (primitive->density * vector_norm2(primitive->velocity) / 2));
}

static void source(void *source_, const double *property, Vector center, double time,
                   const void *context)
{
    (void)context;
    NavierStokesConserved *source = source_;
    double gamma = property[NAVIERSTOKES_HEAT_CAPACITY_RATIO];
    double viscosity = property[NAVIERSTOKES_DYNAMIC_VISCOSITY];
    double prandtl = property[NAVIERSTOKES_PRANDTL];
    double phase = (b * (center.x + center.y + center.z)) - (c * time);
    double sin_phase = sin(phase);
    double cos_phase = cos(phase);
    source->density = a * ((3 * b) - c) * cos_phase;
    source->momentum.x = source->momentum.y = source->momentum.z =
        a * ((b * ((gamma - 1) * ((4 * a * sin_phase) + 5))) + (6 * b) - (2 * c)) * cos_phase / 2;
    source->energy =
        0.5 * a *
        ((2 * a * prandtl * ((3 * b * gamma) - c) * sin(2 * phase)) +
         (6 * b * b * gamma * viscosity * sin_phase) +
         (((15 * b * gamma * prandtl) + (9 * b * prandtl) - (8 * c * prandtl)) * cos_phase)) /
        prandtl;
}

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);

    Vector min_coord = {-1, -1, -1};
    Vector max_coord = {1, 1, 1};
    Triple num_cells = {50, 50, 50};
    Triple periodic = {1, 1, 1};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, periodic);
    mesh_generate(mesh);
    mesh_summary(mesh);

    Equations *eqns = navierstokes_create(mesh);
    equations_create_reference(eqns, exact);
    equations_set_limiter(eqns, "venkatakrishnan", 1);
    equations_set_initial_condition(eqns, "domain", exact, 0);
    equations_set_source(eqns, source);
    equations_set_property(eqns, "dynamic viscosity", 0.1);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 1);
    simulation_set_out_time(sim, 0.1);
    simulation_summary(sim);

    double time = simulation_run(sim);
    simulation_error(sim, time, 0);

    simulation_destroy(sim);
    equations_destroy(eqns);
    mesh_destroy(mesh);

    teal_deinit();
}
