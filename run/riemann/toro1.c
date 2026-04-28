#include <stdio.h>
#include <stdlib.h>

#include "euler.h"

static void exact(void *primitive_, const double *property, Vector center, double time,
                  const void *context)
{
    EulerPrimitive *primitive = primitive_;
    EulerPrimitive left = {.density = 1, .velocity = {.x = 0.75}, .pressure = 1};
    EulerPrimitive right = {.density = 0.125, .pressure = 0.1};
    double gamma = property[EULER_HEAT_CAPACITY_RATIO];
    double loc = (center.x - 0.3) / time;
    *primitive = euler_riemann(left, right, gamma, loc);
}

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);

    char *flux = (argc > 1) ? argv[1] : "hllc";
    int space_order = (argc > 2) ? (int)strtol(argv[2], 0, 10) : 2;

    Vector min_coord = {.x = 0};
    Vector max_coord = {.x = 1};
    Triple num_cells = {.x = 100};
    Triple periodic = {.x = 0};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, periodic);
    mesh_generate(mesh);
    mesh_summary(mesh);

    Equations *eqns = euler_create(mesh);
    equations_create_reference(eqns, exact);
    equations_set_space_order(eqns, space_order);
    equations_set_convective_flux(eqns, flux);
    equations_set_boundary_condition(eqns, "left", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "right", "supersonic outflow", 0, 0);
    equations_set_initial_condition(eqns, "domain", exact, 0);
    equations_summary(eqns);

    String prefix;
    sprintf(prefix, "%s_%s_%d", argv[0], flux, space_order);

    Simulation *sim = simulation_create(eqns, prefix);
    simulation_set_max_time(sim, 0.2);
    simulation_summary(sim);

    double time = simulation_run(sim);
    simulation_error(sim, time, 0);

    simulation_destroy(sim);
    equations_destroy(eqns);
    mesh_destroy(mesh);

    teal_deinit();
}
