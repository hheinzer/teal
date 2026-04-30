#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "navierstokes.h"
#include "sync.h"

static double mach = 0.1, reynolds = 1600;
static void initial(void *variable_, const double *property, Vector center, double time,
                    const void *context);

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);
    int num = (argc > 1) ? (int)strtol(argv[1], 0, 10) : 64;

    Vector min_coord = {-M_PI, -M_PI, -M_PI};
    Vector max_coord = {M_PI, M_PI, M_PI};
    Triple num_cells = {num, num, num};
    Triple periodic = {1, 1, 1};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, periodic);
    mesh_generate(mesh);
    mesh_summary(mesh);

    Equations *eqns = navierstokes_create(mesh);
    equations_set_limiter(eqns, "venkatakrishnan", 1);
    equations_set_initial_condition(eqns, "domain", initial, 0);
    equations_set_property(eqns, "dynamic viscosity", 1 / reynolds);
    equations_summary(eqns);

    char prefix[128];
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

static void initial(void *variable_, const double *property, Vector center, double time,
                    const void *context)
{
    NavierStokesPrimitive *variable = variable_;
    double gamma = property[NAVIERSTOKES_HEAT_CAPACITY_RATIO];
    double pressure0 = 1 / (gamma * mach * mach);
    variable->velocity.x = sin(center.x) * cos(center.y) * cos(center.z);
    variable->velocity.y = (-cos(center.x)) * sin(center.y) * cos(center.z);
    variable->pressure =
        pressure0 + ((cos(2 * center.x) + cos(2 * center.y)) * (cos(2 * center.z) + 2) / 16);
    variable->density = variable->pressure / pressure0;
}
