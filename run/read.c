#include "mesh.h"
#include "teal.h"

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);

    Mesh *mesh = mesh_read("run/mesh.msh");

    mesh_generate(mesh);
    mesh_summary(mesh);

    mesh_write(mesh, argv[0]);

    mesh_free(mesh);
    teal_deinit();
}
