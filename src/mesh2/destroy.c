#include <assert.h>

#include "private.h"
#include "teal2.h"

void mesh2_destroy(Mesh2 *mesh)
{
    assert(mesh);

    teal2_free(mesh->nodes.global);
    teal2_free(mesh->nodes.coord);

    teal2_free(mesh->cells.node.off);
    teal2_free(mesh->cells.node.idx);
    teal2_free(mesh->cells.cell.off);
    teal2_free(mesh->cells.cell.idx);
    teal2_free(mesh->cells.volume);
    teal2_free(mesh->cells.center);
    teal2_free(mesh->cells.projection);
    teal2_free(mesh->cells.offset);

    teal2_free(mesh->faces.node.off);
    teal2_free(mesh->faces.node.idx);
    teal2_free(mesh->faces.cell_idx);
    teal2_free(mesh->faces.area);
    teal2_free(mesh->faces.center);
    teal2_free(mesh->faces.basis);
    teal2_free(mesh->faces.weight);
    teal2_free(mesh->faces.offset);
    teal2_free(mesh->faces.correction);

    teal2_free(mesh->entities.name);
    teal2_free(mesh->entities.cell_off);
    teal2_free(mesh->entities.face_off);
    teal2_free(mesh->entities.rotation);
    teal2_free(mesh->entities.translation);

    teal2_free(mesh->neighbors.tag);
    teal2_free(mesh->neighbors.rank);
    teal2_free(mesh->neighbors.recv_off);
    teal2_free(mesh->neighbors.send.off);
    teal2_free(mesh->neighbors.send.idx);

    teal2_free(mesh);
}
