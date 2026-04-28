#include <math.h>

#include "equations.h"
#include "mesh.h"
#include "teal.h"

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
    teal_init(&argc, &argv);

    Vector min_coord = {0, 0, 0};
    Vector max_coord = {1, 1, 1};
    Triple num_cells = {32, 32, 32};
    Triple periodic = {1, 1, 1};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, periodic);
    mesh_generate(mesh);
    mesh_validate(mesh);
    mesh_summary(mesh);

    Equations *eqns = equations_create(mesh, "test", timestep, 0, 0, 0);
    equations_create_primitive(eqns, (const char *[]){"velocity"}, (int[]){3}, 0, 1);
    equations_set_initial_condition(eqns, "domain", compute, 0);
    equations_summary(eqns);

    mesh_write(mesh, argv[0]);
    equations_write(eqns, argv[0], 0, 0);

    mesh_destroy(mesh);
    equations_destroy(eqns);

    teal_deinit();
}
