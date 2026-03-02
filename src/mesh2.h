#pragma once

#include "utils2.h"
#include "vector2.h"

typedef struct mesh Mesh2;

// Create a distributed Cartesian mesh with optional per-axis periodicity.
Mesh2 *mesh2_create(Vector min_coord, Vector max_coord, Triple num_cells, Triple periodic);

// Read and partition a mesh from disk.
Mesh2 *mesh2_read(const char *fname);

// Modify all node coordinates.
void mesh2_modify(Mesh2 *mesh, void (*modify)(Vector *coord));

// Split an entity by a plane through `root` with normal `normal`.
void mesh2_split(Mesh2 *mesh, const char *entity, Vector root, Vector normal);

// Generate derived connectivity/geometry fields.
void mesh2_generate(Mesh2 *mesh);

// Validate mesh consistency.
void mesh2_validate(const Mesh2 *mesh);

// Print a global mesh summary.
void mesh2_summary(const Mesh2 *mesh);

// Write mesh and derived data to disk.
void mesh2_write(const Mesh2 *mesh, const char *name);

// Release all mesh resources.
void mesh2_destroy(Mesh2 *mesh);
