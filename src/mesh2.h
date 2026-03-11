#pragma once

#include "matrix2.h"
#include "utils2.h"
#include "vector2.h"

enum MeshLimits {
    MAX_CELL_NODES = 8,
    MAX_CELL_FACES = 6,
    MAX_FACE_NODES = 4,
};

typedef struct {
    int num;
    int num_inner;
    long *global;
    Vector *coord;
} MeshNodes;

typedef struct {
    int num;
    int num_inner;
    int off_boundary;
    int off_periodic;
    Graph node;
    Graph cell;
    double *volume, sum_volume;
    Vector *center;
    Vector *projection;
    Vector *offset;
} MeshCells;

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
} MeshFaces;

typedef struct {
    int num;
    int num_inner;
    int off_boundary;
    String *name;
    int *cell_off;
    int *face_off;
    Matrix *rotation;
    Vector *translation;
} MeshEntities;

typedef struct {
    int num;
    int (*tag)[2];
    int *rank;
    int *recv_off;
    Graph send;
} MeshNeighbors;

typedef struct {
    MeshNodes nodes;
    MeshCells cells;
    MeshFaces faces;
    MeshEntities entities;
    MeshNeighbors neighbors;
} Mesh;

// Create a distributed Cartesian mesh with optional per-axis periodicity.
Mesh *mesh2_create(Vector min_coord, Vector max_coord, Triple num_cells, Triple periodic);

// Read and partition a mesh from disk.
Mesh *mesh2_read(const char *fname);

// Modify all node coordinates.
void mesh2_modify(Mesh *mesh, void (*modify)(Vector *coord));

// Split an entity by a plane through `root` with normal `normal`.
void mesh2_split(Mesh *mesh, const char *entity, Vector root, Vector normal);

// Generate derived connectivity/geometry fields.
void mesh2_generate(Mesh *mesh);

// Validate mesh consistency.
void mesh2_validate(const Mesh *mesh);

// Print a global mesh summary.
void mesh2_summary(const Mesh *mesh);

// Write mesh and derived data to disk.
void mesh2_write(const Mesh *mesh, const char *name);

// Release all mesh resources.
void mesh2_destroy(Mesh *mesh);
