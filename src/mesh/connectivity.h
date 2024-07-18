#ifndef CONNECTIVITY_H
#define CONNECTIVITY_H

#include "mesh.h"

void connectivity_cells(Mesh *mesh);

void connectivity_faces(Mesh *mesh);

long connectivity_nodes(const Mesh *mesh, long *node, long j0, long j1);

#endif
