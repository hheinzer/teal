#include <math.h>
#include <stdio.h>

#include "mesh.h"

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);

    Mesh *mesh = mesh_read("run/test/mesh.h5");

    for (long i = 0; i < mesh->nodes.num; i++) {
        vector *coord = &mesh->nodes.coord[i];
        scalar a = 4.5, b = 3.5, c = 1.0 / 6;
        coord->y += c * (3 - coord->y) * (1 + tanh(b * (fabs(coord->x - a) - b)));
    }

    vector root = {3, 0, 0};
    vector normal = {1, 0, 0};
    mesh_split(mesh, "bottom", root, normal);

    mesh_commit(mesh);

    mesh_check(mesh);
    mesh_summary(mesh);

    char fname[128];
    sprintf(fname, "%s_mesh.h5", argv[0]);
    mesh_write(mesh, fname);

    teal_finalize();
}
