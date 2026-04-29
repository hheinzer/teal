#include <math.h>

#include "euler.h"

static double wrap(double x, double xmin, double xmax)
{
    return (x >= 0 ? xmin : xmax) + fmod(x, xmax - xmin);
}

static void exact(void *primitive_, const double *property, Vector center, double time,
                  const void *context)
{
    (void)context;
    EulerPrimitive *primitive = primitive_;
    double gamma = property[EULER_HEAT_CAPACITY_RATIO];
    Vector position = {5, 5, 0};
    Vector velocity = {1, 1, 0};
    double vortex_strength = 5;
    double dx = wrap(center.x - (velocity.x * time), 0, 10) - position.x;
    double dy = wrap(center.y - (velocity.y * time), 0, 10) - position.y;
    double r2 = (dx * dx) + (dy * dy);
    primitive->density =
        pow(1 - ((gamma - 1) * sq(vortex_strength) / (8 * gamma * sq(M_PI)) * exp(1 - r2)),
            1 / (gamma - 1));
    primitive->velocity.x = velocity.x - (vortex_strength / (2 * M_PI) * exp((1 - r2) / 2) * dy);
    primitive->velocity.y = velocity.y + (vortex_strength / (2 * M_PI) * exp((1 - r2) / 2) * dx);
    primitive->velocity.z = 0;
    primitive->pressure = pow(primitive->density, gamma);
}

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);

    Vector min_coord = {0, 0, 0};
    Vector max_coord = {10, 10, 0};
    Triple num_cells = {256, 256, 0};
    Triple periodic = {1, 1, 0};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, periodic);
    mesh_generate(mesh);
    mesh_summary(mesh);

    Equations *eqns = euler_create(mesh);
    equations_create_reference(eqns, exact);
    equations_set_initial_condition(eqns, "domain", exact, 0);
    equations_summary(eqns);

    Simulation *sim = simulation_create(eqns, argv[0]);
    simulation_set_max_time(sim, 10);
    simulation_set_out_time(sim, 1);
    simulation_summary(sim);

    double time = simulation_run(sim);
    simulation_error(sim, time, 0);

    simulation_destroy(sim);
    equations_destroy(eqns);
    mesh_destroy(mesh);

    teal_deinit();
}
