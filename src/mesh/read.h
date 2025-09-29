#pragma once

#include "mesh.h"

void mesh_read_hdf5(Mesh *mesh, const char *fname);

void mesh_read_gmsh(Mesh *mesh, const char *fname);
