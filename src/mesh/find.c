#include "find.h"

#include <stdlib.h>
#include <string.h>

long mesh_find_entity(const Mesh *mesh, const char *entity)
{
    for (long i = 0; i < mesh->n_entities; ++i)
        if (!strcmp(mesh->entity.name[i], entity)) return i;
    abort();
}
