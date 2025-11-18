#pragma once

#include "mesh.h"

/* Reorder all nodes by map; remap cell->node if present. */
void mesh_reorder_nodes(MeshNodes *nodes, MeshCells *cells, const int *map);

/* Reorder cells in [beg,end) by map; remap cell->cell and faces->cell if present. */
void mesh_reorder_cells(MeshCells *cells, MeshFaces *faces, int beg, int end, const int *map);

/* Reorder all faces by key. */
void mesh_reorder_faces(MeshFaces *faces, const int *key);
