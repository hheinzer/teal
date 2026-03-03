#include <assert.h>

#include "mesh2.h"

void mesh2_modify(Mesh *mesh, void (*modify)(Vector *))
{
    assert(mesh && modify);
    for (int i = 0; i < mesh->nodes.num; i++) {
        modify(&mesh->nodes.coord[i]);
    }
}
