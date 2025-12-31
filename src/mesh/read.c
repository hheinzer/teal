#include <assert.h>
#include <string.h>

#include "mesh.h"
#include "read_gmsh.h"
#include "teal.h"

// Read a mesh file based on its extension.
static void read_file(Mesh *mesh, const char *fname)
{
    char *ext = strrchr(fname, '.');
    if (!ext) {
        teal_error("missing file extension (%s)", fname);
    }
    else if (!strcmp(ext, ".msh")) {
        read_gmsh(mesh, fname);
    }
    else {
        teal_error("invalid file extension (%s)", ext);
    }
}

Mesh *mesh_read(const char *fname)
{
    assert(fname);
    Mesh *mesh = teal_alloc(1, sizeof(*mesh));

    read_file(mesh, fname);

    return mesh;
}
