#pragma once

#include "mesh.h"

long connectivity_nodes(const Mesh *mesh, long *node, long cell_a, long cell_b);

void connectivity_cells(Mesh *mesh);

void connectivity_faces(Mesh *mesh);
