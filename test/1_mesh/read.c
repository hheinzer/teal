#include <math.h>

#include "mesh.h"
#include "teal.h"

static void modify(Vector *coord)
{
    double a = 4.5, b = 3.5, c = 6.0;
    coord->y += (3 - coord->y) * (1 + tanh(b * (fabs(coord->x - a) - b))) / c;
}

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);

    Mesh *mesh = mesh_read("test/1_mesh/mesh.msh");

    mesh_modify(mesh, modify);

    Vector root = {3, 0, 0};
    Vector normal = {1, 0, 0};
    mesh_split(mesh, "bottom", root, normal);

    mesh_generate(mesh);
    mesh_validate(mesh);
    mesh_summary(mesh);

    mesh_write(mesh, argv[0]);

    mesh_destroy(mesh);

    teal_deinit();
}
