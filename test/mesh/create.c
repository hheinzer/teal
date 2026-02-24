#include <math.h>

#include "mesh2.h"
#include "teal2.h"

static void modify(Vector *coord)
{
    double a = 4.5, b = 3.5, c = 6.0;
    coord->y += (3 - coord->y) * (1 + tanh(b * (fabs(coord->x - a) - b))) / c;
}

int main(int argc, char **argv)
{
    teal2_init(&argc, &argv);

    Vector min_coord = {0, 0, 0};
    Vector max_coord = {9, 3, 3};
    Triple num_cells = {30, 10, 10};
    Triple periodic = {1, 0, 1};
    Mesh2 *mesh = mesh2_create(min_coord, max_coord, num_cells, periodic);

    mesh2_modify(mesh, modify);

    Vector root = {3, 0, 0};
    Vector normal = {1, 0, 0};
    mesh2_split(mesh, "bottom", root, normal);

    mesh2_generate(mesh);
    mesh2_validate(mesh);

    mesh2_summary(mesh);
    mesh2_write(mesh, argv[0]);

    mesh2_destroy(mesh);

    teal2_deinit();
}
