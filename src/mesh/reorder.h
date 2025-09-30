#pragma once

#include "mesh.h"

/* Reorder all nodes by map; remap cell->node if present. */
void mesh_reorder_nodes(MeshNodes *nodes, MeshCells *cells, const long *map);

/* Reorder cells in [beg,end) by map; remap cell->cell and faces->cell if present. */
void mesh_reorder_cells(MeshCells *cells, MeshFaces *faces, long beg, long end, const long *map);

/* Reorder all faces by key. */
void mesh_reorder_faces(MeshFaces *faces, const long *key);
