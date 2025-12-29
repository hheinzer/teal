#include <assert.h>

#include "mesh.h"
#include "teal.h"

void mesh_free(Mesh *mesh)
{
    assert(mesh);
    teal_free(mesh);
}
