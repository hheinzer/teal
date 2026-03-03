#pragma once

#include "matrix2.h"
#include "utils2.h"
#include "vector2.h"

enum Mesh2Limits {
    MAX_CELL_NODES = 8,
    MAX_CELL_FACES = 6,
    MAX_FACE_NODES = 4,
};

typedef struct {
    int num;
    int num_inner;
    long *global;
    Vector *coord;
} Mesh2Nodes;

typedef struct {
    int num;
    int num_inner;
    int off_boundary;
    int off_periodic;
    Graph node;
    Graph cell;
    double *volume;
    Vector *center;
    Vector *projection;
    Vector *offset;
} Mesh2Cells;

typedef struct {
    int num;
    int num_inner;
    int off_boundary;
    Graph node;
    Pair *cell_idx;
    double *area;
    Vector *center;
    Matrix *basis;
    VectorPair *weight;
    VectorPair *offset;
    VectorNorm *correction;
} Mesh2Faces;

typedef struct {
    int num;
    int num_inner;
    int off_boundary;
    String *name;
    int *cell_off;
    int *face_off;
    Matrix *rotation;
    Vector *translation;
} Mesh2Entities;

typedef struct {
    int num;
    int (*tag)[2];
    int *rank;
    int *recv_off;
    Graph send;
} Mesh2Neighbors;

typedef struct {
    Mesh2Nodes nodes;
    Mesh2Cells cells;
    Mesh2Faces faces;
    Mesh2Entities entities;
    Mesh2Neighbors neighbors;
} Mesh2;

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
