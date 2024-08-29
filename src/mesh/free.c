#include "mesh.h"
#include "teal/memory.h"

void mesh_free(Mesh **mesh)
{
    memory_free(&(*mesh)->node.global);
    memory_free(&(*mesh)->node.coord);

    memory_free(&(*mesh)->cell.i_node);
    memory_free(&(*mesh)->cell.node);
    memory_free(&(*mesh)->cell.i_cell);
    memory_free(&(*mesh)->cell.cell);
    memory_free(&(*mesh)->cell.volume);
    memory_free(&(*mesh)->cell.center);
    memory_free(&(*mesh)->cell.projection);
    memory_free(&(*mesh)->cell.to_cell);

    memory_free(&(*mesh)->face.i_node);
    memory_free(&(*mesh)->face.node);
    memory_free(&(*mesh)->face.cell);
    memory_free(&(*mesh)->face.area);
    memory_free(&(*mesh)->face.center);
    memory_free(&(*mesh)->face.basis);
    memory_free(&(*mesh)->face.to_cell);
    memory_free(&(*mesh)->face.weight);
    memory_free(&(*mesh)->face.correction);

    memory_free(&(*mesh)->entity.name);
    memory_free(&(*mesh)->entity.j_cell);
    memory_free(&(*mesh)->entity.j_face);
    memory_free(&(*mesh)->entity.offset);

    memory_free(&(*mesh)->sync.j_recv);
    memory_free(&(*mesh)->sync.i_send);
    memory_free(&(*mesh)->sync.send);

    memory_free(mesh);
}
