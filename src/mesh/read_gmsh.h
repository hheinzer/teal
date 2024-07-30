#pragma once

#include "mesh.h"

/* Create a 'mesh' from a gmsh ".geo" or ".msh" file 'fname'. */
void read_gmsh(Mesh *mesh, const char *fname);
