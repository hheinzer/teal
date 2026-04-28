#pragma once

#include "mesh.h"

void mesh_reorder_nodes(Mesh *mesh, const int *map, int beg, int end);

void mesh_reorder_cells(Mesh *mesh, const int *map, int beg, int end);
