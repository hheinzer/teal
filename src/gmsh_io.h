#ifndef GMSH_IO_H
#define GMSH_IO_H

#include "mesh.h"

void gmsh_read(Mesh *mesh, const char *fname);

void gmsh_export(Mesh *mesh);

#endif
