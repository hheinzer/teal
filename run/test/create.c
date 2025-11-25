#include <math.h>

#include "mesh.h"

int main(int argc, char **argv)
{
    teal_initialize(&argc, &argv);

    vector min_coord = {0, 0, 0};
    vector max_coord = {9, 3, 1};
    tuple num_cells = {300, 100, 33};
    bool periodic[3] = {true, false, true};
    Mesh *mesh = mesh_create(min_coord, max_coord, num_cells, periodic, 3);

    for (long i = 0; i < mesh->nodes.num; i++) {
        vector *coord = &mesh->nodes.coord[i];
        scalar a = 4.5, b = 3.5, c = 1.0 / 6;
        coord->y += c * (3 - coord->y) * (1 + tanh(b * (fabs(coord->x - a) - b)));
    }

    vector root = {3, 0, 0};
    vector normal = {1, 0, 0};
    mesh_split(mesh, "bottom", root, normal);

    mesh_generate(mesh);

    mesh_check(mesh);
    mesh_summary(mesh);

    mesh_write(mesh, argv[0]);

    teal_finalize();
}
