#include "mesh.h"

void mesh_modify(Mesh *mesh, Modify *modify)
{
    for (long i = 0; i < mesh->n_nodes; ++i) modify(mesh->node.coord[i]);
}
