#pragma once

#include "mesh2.h"

void mesh2_reorder_nodes(Mesh *mesh, const int *map, int beg, int end);

void mesh2_reorder_cells(Mesh *mesh, const int *map, int beg, int end);
