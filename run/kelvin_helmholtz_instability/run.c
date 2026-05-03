#include <math.h>

#include "euler.h"

static void initial(void *primitive_, const double *property, Vector center, double time,
                    const void *context)
{
    EulerPrimitive *primitive = primitive_;
    double w0 = 0.1, sigma2 = sq(0.05 / 2);
    if (0.25 < center.y && center.y < 0.75) {
        primitive->density = 2;
        primitive->velocity.x = 0.5;
    }
    else {
        primitive->density = 1;
        primitive->velocity.x = -0.5;
    }
    primitive->velocity.y =
        w0 * sin(4 * M_PI * center.x) *
        (exp(-sq(center.y - 0.25) / (2 * sigma2)) + exp(-sq(center.y - 0.75) / (2 * sigma2)));
    primitive->pressure = 2.5;
}

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);

    Vector min_coord = {.x = 0, .y = 0};
    Vector max_coord = {.x = 1, .y = 1};
    Triple num_cells = {.x = 256, .y = 256};
    Triple periodic = {.x = 1, .y = 1};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, periodic);
    mesh_generate(mesh);
    mesh_summary(mesh);

    Equations *eqns = euler_create(mesh);
    equations_set_limiter(eqns, "venkatakrishnan", 1);
    equations_set_initial_condition(eqns, "domain", initial, 0);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 1);
    simulation_set_out_time(sim, 0.1);
    simulation_summary(sim);

    simulation_run(sim);

    simulation_destroy(sim);
    equations_destroy(eqns);
    mesh_destroy(mesh);

    teal_deinit();
}
