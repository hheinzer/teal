#include <stdlib.h>
#include <string.h>

#include "connectivity.h"
#include "mesh.h"
#include "read_gmsh.h"
#include "teal/memory.h"
#include "teal/sync.h"

Mesh *mesh_read(const char *fname)
{
    memory_sum_setzero();

    Mesh *mesh = memory_calloc(1, sizeof(*mesh));
    if (sync.rank != 0) return mesh;

    const char *ext = strrchr(fname, '.');
    if (!strcmp(ext, ".geo") || !strcmp(ext, ".msh"))
        read_gmsh(mesh, fname);
    else
        abort();

    connectivity_cells(mesh);

    return mesh;
}
