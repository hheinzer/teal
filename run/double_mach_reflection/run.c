#include <math.h>

#include "euler.h"

static const EulerPrimitive upstream = {.density = 1.4, .pressure = 1};
static const EulerPrimitive downstream = {
    .density = 8, .velocity = {.x = 7.1447096, .y = -4.125}, .pressure = 116.5};
static const double alpha = 30 * M_PI / 180, shift = 1.0 / 6, Wx = 8.660254, Wy = -5;

static void initial(void *primitive_, const double *property, Vector center, double time,
                    const void *context)
{
    EulerPrimitive *primitive = primitive_;
    double angle = atan2(center.x - shift - (Wx * time), center.y - (Wy * time));
    *primitive = (angle >= alpha) ? upstream : downstream;
}

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);

    Vector min_coord = {.x = 0, .y = 0};
    Vector max_coord = {.x = 3.25, .y = 1};
    Triple num_cells = {.x = 1040, .y = 320};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, (Triple){0});

    Vector root = {.x = shift};
    Vector normal = {.x = 1};
    mesh_split(mesh, "bottom", root, normal);

    mesh_generate(mesh);
    mesh_summary(mesh);

    Equations *eqns = euler_create(mesh);
    equations_set_boundary_condition(eqns, "left", "supersonic inflow", 0, &downstream);
    equations_set_boundary_condition(eqns, "right", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "bottom-a", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "bottom-b", "slipwall", 0, 0);
    equations_set_boundary_condition(eqns, "top", "custom", initial, 0);
    equations_set_initial_condition(eqns, "domain", initial, 0);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 0.2);
    simulation_set_out_time(sim, 0.02);
    simulation_summary(sim);

    simulation_run(sim);

    simulation_destroy(sim);
    equations_destroy(eqns);
    mesh_destroy(mesh);

    teal_deinit();
}
