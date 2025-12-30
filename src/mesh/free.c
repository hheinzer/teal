#include "mesh.h"
#include "teal.h"

void mesh_free(Mesh *mesh)
{
    if (!mesh) {
        return;
    }

    teal_free(mesh->nodes.global);
    teal_free(mesh->nodes.coord);

    teal_free(mesh->cells.node.off);
    teal_free(mesh->cells.node.idx);
    teal_free(mesh->cells.cell.off);
    teal_free(mesh->cells.cell.idx);
    teal_free(mesh->cells.volume);
    teal_free(mesh->cells.center);
    teal_free(mesh->cells.projection);
    teal_free(mesh->cells.offset);

    teal_free(mesh->faces.node.off);
    teal_free(mesh->faces.node.idx);
    teal_free(mesh->faces.cell);
    teal_free(mesh->faces.area);
    teal_free(mesh->faces.center);
    teal_free(mesh->faces.basis);
    teal_free(mesh->faces.weight);
    teal_free(mesh->faces.offset);
    teal_free(mesh->faces.correction);

    teal_free(mesh->entities.name);
    teal_free(mesh->entities.cell_off);
    teal_free(mesh->entities.face_off);

    teal_free(mesh->periodics.cell_off);
    teal_free(mesh->periodics.face_off);
    teal_free(mesh->periodics.translation);

    teal_free(mesh->neighbors.rank);
    teal_free(mesh->neighbors.recv_off);
    teal_free(mesh->neighbors.send.off);
    teal_free(mesh->neighbors.send.idx);

    teal_free(mesh);
}
