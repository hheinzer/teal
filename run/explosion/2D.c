#include "euler.h"

static void initial(void *primitive_, const double *property, Vector center, double time,
                    const void *context)
{
    EulerPrimitive *primitive = primitive_;
    EulerPrimitive inner = {.density = 1, .pressure = 1};
    EulerPrimitive outer = {.density = 0.125, .pressure = 0.1};
    *primitive = (vector_norm(center) <= 0.4) ? inner : outer;
}

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);

    Vector min_coord = {.x = -1, .y = -1, .z = -1};
    Vector max_coord = {.x = 1, .y = 1, .z = 1};
    Triple num_cells = {.x = 100, .y = 100};
    Triple periodic = {.x = 0};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, periodic);
    mesh_generate(mesh);
    mesh_summary(mesh);

    Equations *eqns = euler_create(mesh);
    equations_set_boundary_condition(eqns, "left", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "right", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "bottom", "supersonic outflow", 0, 0);
    equations_set_boundary_condition(eqns, "top", "supersonic outflow", 0, 0);
    equations_set_initial_condition(eqns, "domain", initial, 0);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 0.25);
    simulation_summary(sim);

    simulation_run(sim);

    simulation_destroy(sim);
    equations_destroy(eqns);
    mesh_destroy(mesh);

    teal_deinit();
}
