#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "navierstokes.h"
#include "sync.h"

static const double mach = 0.1, reynolds = 1600;

static void initial(void *primitive_, const double *property, Vector center, double time,
                    const void *context)
{
    NavierStokesPrimitive *primitive = primitive_;
    double gamma = property[NAVIERSTOKES_HEAT_CAPACITY_RATIO];
    double pressure0 = 1 / (gamma * sq(mach));
    primitive->velocity.x = sin(center.x) * cos(center.y) * cos(center.z);
    primitive->velocity.y = -cos(center.x) * sin(center.y) * cos(center.z);
    primitive->pressure =
        pressure0 + ((cos(2 * center.x) + cos(2 * center.y)) * (cos(2 * center.z) + 2) / 16);
    primitive->density = primitive->pressure / pressure0;
}

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);

    int num = (argc > 1) ? (int)strtol(argv[1], 0, 10) : 64;

    Vector min_coord = {.x = -M_PI, .y = -M_PI, .z = -M_PI};
    Vector max_coord = {.x = M_PI, .y = M_PI, .z = M_PI};
    Triple num_cells = {.x = num, .y = num, .z = num};
    Triple periodic = {.x = 1, .y = 1, .z = 1};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, periodic);
    mesh_generate(mesh);
    mesh_summary(mesh);

    Equations *eqns = navierstokes_create(mesh);
    equations_set_limiter(eqns, "venkatakrishnan", 1);
    equations_set_initial_condition(eqns, "domain", initial, 0);
    equations_set_property(eqns, "dynamic viscosity", 1 / reynolds);
    equations_summary(eqns);

    String prefix;
    sprintf(prefix, "%s_%04d_%04d", argv[0], num, sync.size);

    Simulation *sim = simulation_create(eqns, prefix);
    simulation_set_max_time(sim, 1);
    simulation_set_out_time(sim, 0.1);
    simulation_summary(sim);

    simulation_run(sim);

    simulation_destroy(sim);
    equations_destroy(eqns);
    mesh_destroy(mesh);

    teal_deinit();
}
