#include <stdlib.h>

#include "mesh.h"

void mesh_free(Mesh *mesh)
{
    free(mesh->node.idx);
    free(mesh->node.coord);

    free(mesh->cell.i_node);
    free(mesh->cell.node);
    free(mesh->cell.i_cell);
    free(mesh->cell.cell);
    free(mesh->cell.volume);
    free(mesh->cell.center);
    free(mesh->cell.projection);
    free(mesh->cell.reconstruction);

    free(mesh->face.i_node);
    free(mesh->face.node);
    free(mesh->face.cell);
    free(mesh->face.area);
    free(mesh->face.center);
    free(mesh->face.normal);
    free(mesh->face.gradient_weight);
    free(mesh->face.reconstruction);
    free(mesh->face.gradient_correction);

    free(mesh->entity.name);
    free(mesh->entity.j_cell);
    free(mesh->entity.j_face);
    free(mesh->entity.offset);

    free(mesh->sync.j_recv);
    free(mesh->sync.i_send);
    free(mesh->sync.send);

    *mesh = (Mesh){0};
}
