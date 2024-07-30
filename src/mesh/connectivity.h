#pragma once

#include "mesh.h"

/* Compute the cell to cell map of the 'mesh'. */
void connectivity_cells(Mesh *mesh);

/* Compute the face to node map, the face to cell map, and the entity to face map of the 'mesh'. */
void connectivity_faces(Mesh *mesh);

/* Compute the shared 'node' indices of the two cells 'j0' and 'j1' of the 'mesh'. WARNING: This
 * function only works for unique node indices. */
long connectivity_nodes(const Mesh *mesh, long *node, long j0, long j1);
