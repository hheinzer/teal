#include <assert.h>

#include "private.h"

void mesh2_modify(Mesh2 *mesh, void (*modify)(Vector *))
{
    assert(mesh && !mesh->generated && modify);
    for (int i = 0; i < mesh->nodes.num; i++) {
        modify(&mesh->nodes.coord[i]);
    }
}
