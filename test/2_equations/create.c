#include <math.h>

#include "equations2.h"
#include "mesh2.h"
#include "teal2.h"

static double timestep(const void *primitive, const double *property, double volume,
                       Vector projection)
{
    return 1;
}

static void compute(void *primitive_, const double *property, Vector center, double time,
                    const void *context)
{
    double fac = 2 * M_PI;
    Vector *velocity = primitive_;
    velocity->x = sin(fac * center.x) * cos(fac * center.y) * cos(fac * center.z);
    velocity->y = cos(fac * center.x) * sin(fac * center.y) * cos(fac * center.z);
    velocity->z = 0;
}

int main(int argc, char **argv)
{
    teal2_init(&argc, &argv);

    Vector min_coord = {0, 0, 0};
    Vector max_coord = {1, 1, 1};
    Triple num_cells = {32, 32, 32};
    Triple periodic = {1, 1, 1};
    Mesh *mesh = mesh2_create(min_coord, max_coord, num_cells, periodic);
    mesh2_generate(mesh);
    mesh2_validate(mesh);
    mesh2_summary(mesh);

    Equations *eqns = equations2_create(mesh, "test", timestep, 0, 0, 0);
    equations2_create_primitive(eqns, (const char *[]){"velocity"}, (int[]){3}, 0, 1);
    equations2_set_initial_condition(eqns, "domain", compute, 0);
    equations2_summary(eqns);

    mesh2_write(mesh, argv[0]);
    equations2_write(eqns, argv[0], 0, 0);

    mesh2_destroy(mesh);
    equations2_destroy(eqns);

    teal2_deinit();
}
