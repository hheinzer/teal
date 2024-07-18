#include <string.h>

#include "connectivity.h"
#include "core/utils.h"
#include "mesh.h"
#include "read_gmsh.h"
#include "teal.h"

Mesh mesh_read(const char *fname)
{
    Mesh mesh = {0};
    if (teal.rank != 0) return mesh;

    const char *ext = strrchr(fname, '.');
    if (!strcmp(ext, ".geo") || !strcmp(ext, ".msh"))
        read_gmsh(&mesh, fname);
    else
        error("unsupported mesh format '%s'", ext);

    connectivity_cells(&mesh);

    return mesh;
}
