#include <math.h>

#include "euler.h"

static void initial(void *primitive_, const double *property, Vector center, double time,
                    const void *context)
{
    EulerPrimitive *primitive = primitive_;
    primitive->pressure = 1;
    if (center.x >= 18 + (2 * cos(2 * M_PI * center.y / 15))) {
        primitive->density = 0.25;
    }
    else {
        primitive->density = 1;
    }
    if (2 <= center.x && center.x <= 6) {
        primitive->density = 4.22;
        primitive->pressure = 4.9;
    }
}

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);

    Vector min_coord = {.x = 0, .y = 0, .z = 0};
    Vector max_coord = {.x = 60, .y = 15, .z = 15};
    Triple num_cells = {.x = 400, .y = 100, .z = 1};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, (Triple){0});
    mesh_generate(mesh);
    mesh_summary(mesh);

    Equations *eqns = euler_create(mesh);
    equations_set_limiter(eqns, "venkatakrishnan", 1);
    equations_set_boundary_condition(eqns, "left", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "right", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "bottom", "symmetry", 0, 0);
    equations_set_boundary_condition(eqns, "top", "symmetry", 0, 0);
    equations_set_initial_condition(eqns, "domain", initial, 0);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 300);
    simulation_set_out_time(sim, 30);
    simulation_summary(sim);

    simulation_run(sim);

    simulation_destroy(sim);
    equations_destroy(eqns);
    mesh_destroy(mesh);

    teal_deinit();
}
